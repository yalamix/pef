from sympy import Float, sympify, Symbol, Piecewise, Eq, simplify, Basic, Rational, Expr, Add, Number
from sympy.assumptions.ask import ask
from sympy.assumptions import Q
from sympy.functions.special.singularity_functions import SingularityFunction
from sympy.printing.latex import LatexPrinter
from sympy import (
    SingularityFunction, symbols, sympify, simplify, Basic
)
from sympy.core.relational import Gt, Ge, Lt, Le, Relational
from typing import List, Union, Iterable, Mapping
from PIL import Image
import plotly.graph_objects as go
import numpy as np
import io
import base64
import re

mapping = {
    'SingularityFunction':  # name that appears in the generated lambda
        lambda x,a,n: np.where(x > a, (x - a)**n, 0)
}

class ThresholdLatexPrinter(LatexPrinter):
    def __init__(self, *args, sci_min=-2, sci_max=4, **kwargs):
        super().__init__(*args, **kwargs)
        self.sci_min = sci_min
        self.sci_max = sci_max

    def _print_Float(self, f: Float) -> str:
        """
        1) If f is within 1e-12 of an integer, print as “175” (no .00).
        2) Else if f can be written as an irreducible fraction with numerator
           and denominator both < 100 (and denominator != 1), print as
           “\frac{num}{den}”.
        3) Else if exponent outside [sci_min, sci_max], use scientific notation
           “m\times10^{e}” with mantissa from "{:.15e}" logic.
        4) Otherwise, round to two decimal places (“%.2f”).
        """
        pyf = float(f)

        # 2) Fraction‐check with max denominator 999 (skip denom=1)
        R = Rational(pyf).limit_denominator(99)
        if R.q != 1 and abs(float(R) - pyf) < 1e-12 \
           and abs(R.p) < 100 and abs(R.q) < 100:
            return rf"\frac{{{R.p}}}{{{R.q}}}"

        # 3) Scientific‐notation check
        sci_str = "{:.3e}".format(pyf)
        mant_str, exp_str = sci_str.split("e")
        e = int(exp_str)
        mant_str = mant_str.rstrip("0").rstrip(".")
        if e > self.sci_max or e < self.sci_min:
            return rf"{mant_str}\times10^{{{e}}}"
        # 1) Integer check via rounding within a tiny epsilon
        elif abs(pyf - round(pyf)) < 1e-12:
            return str(int(round(pyf)))        

        # 4) Otherwise, regular float rounded to 2 decimals
        return f"{pyf:.2f}"

    def _print_SingularityFunction(self, expr: SingularityFunction) -> str:
        """
        Render <x - a>^n as LaTeX.  If `a` is not an atom (e.g. `c+b`), wrap it in parentheses:
          … \left\langle x - (c + b) \right\rangle^{n} …
        Otherwise, if `a` is a plain symbol or number, no parentheses are needed:
          … \left\langle x - a \right\rangle^{n} …
        """
        xx, aa, nn = expr.args

        # 1) Convert x and n in the usual way
        px = self._print(xx)
        pn = self._print(nn)

        # 2) Convert `aa` to LaTeX; then decide if we need parentheses
        pa = self._print(aa)

        # If `aa` is not a single atom (Symbol, Number, etc.), wrap it:
        #   - aa.is_Atom is True for Symbol or Number
        #   - but an Add or Mul or any composite expression has is_Atom=False
        if not aa.is_Atom:
            # Use \bigl(...\bigr) to avoid overly-large delimiters
            pa = r'\bigl(' + pa + r'\bigr)'

        return r'{\left\langle %s - %s \right\rangle}^{%s}' % (px, pa, pn)

    def _print_Symbol(self, expr):
        name = expr.name
        parts = name.split('_')
        base = parts[0]
        subscripts = parts[1:]

        # Mapping of Greek letter names to LaTeX commands
        greek_letters = {
            'alpha': r'\alpha',
            'beta': r'\beta',
            'gamma': r'\gamma',
            'delta': r'\delta',
            'epsilon': r'\epsilon',
            'zeta': r'\zeta',
            'eta': r'\eta',
            'theta': r'\theta',
            'iota': r'\iota',
            'kappa': r'\kappa',
            'lambda': r'\lambda',
            'mu': r'\mu',
            'nu': r'\nu',
            'xi': r'\xi',
            'omicron': r'o',
            'pi': r'\pi',
            'rho': r'\rho',
            'sigma': r'\sigma',
            'tau': r'\tau',
            'upsilon': r'\upsilon',
            'phi': r'\phi',
            'chi': r'\chi',
            'psi': r'\psi',
            'omega': r'\omega',
            'Gamma': r'\Gamma',
            'Delta': r'\Delta',
            'Theta': r'\Theta',
            'Lambda': r'\Lambda',
            'Xi': r'\Xi',
            'Pi': r'\Pi',
            'Sigma': r'\Sigma',
            'Upsilon': r'\Upsilon',
            'Phi': r'\Phi',
            'Psi': r'\Psi',
            'Omega': r'\Omega',
        }

        # Convert base to LaTeX if it's a Greek letter
        base_latex = greek_letters.get(base, base)

        # Construct nested subscripts
        if subscripts:
            subscript = subscripts[-1]
            for s in reversed(subscripts[:-1]):
                subscript = f"{s}_{{{subscript}}}"
            return f"{base_latex}_{{{subscript}}}"
        else:
            return base_latex      

    def _print_Equality(self, expr: Eq) -> str:
        lhs, rhs = expr.lhs, expr.rhs

        # 1) Only box simple Symbol = Number
        if isinstance(lhs, Symbol) and isinstance(rhs, Number):
            lhs_str = self._print(lhs)
            rhs_str = self._print(rhs)
            return r"\boxed{\strut " + lhs_str + r"=" + rhs_str + r"}"

        # 2) Otherwise delegate to the built-in relational printer
        return super()._print_Relational(expr)        

def remove_negative_singularity_terms(expr):
    """
    Removes terms containing SingularityFunction with negative exponents from the expression.
    
    Parameters:
        expr (sympy expression): The input expression.
        
    Returns:
        sympy expression: The expression with negative exponent SingularityFunction terms removed.
    """
    # Expand the expression to separate terms
    expr = expr.expand()
    
    # If the expression is an addition of terms
    if isinstance(expr, Add):
        new_terms = []
        for term in expr.args:
            # Check if the term contains a SingularityFunction with negative exponent
            singularity_funcs = term.atoms(SingularityFunction)
            if any(sf.args[2].is_number and sf.args[2] < 0 for sf in singularity_funcs):
                continue  # Skip this term
            new_terms.append(term)
        return Add(*new_terms)
    else:
        # For non-additive expressions, check directly
        singularity_funcs = expr.atoms(SingularityFunction)
        if any(sf.args[2].is_number and sf.args[2] < 0 for sf in singularity_funcs):
            return 0  # Remove the entire expression
        return expr

def unify_symbols(expr: Union[Basic, Iterable[Basic]],
                  symdict: Mapping[str, Basic]
                 ) -> Union[Basic, list]:
    """
    Replace in `expr` any Symbol whose .name() matches a key in `symdict`
    by the corresponding symdict[name].

    - expr may be a single SymPy expression or an iterable of them.
    - symdict maps symbol-names (str) → the canonical Symbol objects you want.
    Returns the new expression (or list of new expressions), leaving the
    originals untouched.
    """
    def _replace_single(e: Basic) -> Basic:
        # build replacement map: old_symbol → new_symbol
        rep = {
            old: new
            for old in e.free_symbols
            for name, new in symdict.items()
            if old.name == name and old is not new
        }
        return e.xreplace(rep)

    if isinstance(expr, Basic):
        return _replace_single(expr)
    else:
        # assume it's an iterable of expressions
        return [_replace_single(e) for e in expr]

Item = List[Union[str, float, int]]

def add_item(collection: List[Item], new_item: Item) -> bool:
    """
    Add `new_item` to `collection` if there is no existing item in `collection`
    whose first two elements match new_item's first two elements.
    
    Returns True if new_item was added, False otherwise.
    """
    # Compare only the first two elements
    key_new = tuple(new_item[:2])
    for item in collection:
        if tuple(item[:2]) == key_new:
            # Already present (regardless of last element), do not add
            return False
    collection.append(new_item)
    return True

def remove_zero_flags(collection: List[Item]) -> None:
    """
    Remove all items from `collection` whose last element is 0 (or False).
    This mutates the list in place.
    """
    # Iterate backwards so removals don't shift upcoming indices
    for idx in range(len(collection)-1, -1, -1):
        if collection[idx][-1] == 0:
            del collection[idx]


def collapse_singularity_functions(expr: Basic, var: symbols, relations: list) -> Basic:
    """
    Given a SymPy expression `expr` containing SingularityFunction instances,
    collapse every instance of SingularityFunction(f, knot, n) where `f` does not
    contain the variable `var`, based on known `relations` among constants.

    Parameters
    ----------
    expr : sympy.Basic
        The expression in which to find and collapse SingularityFunction instances.
    var : sympy.Symbol
        The “variable” symbol (e.g. x).  If SingularityFunction’s first argument `f`
        has `var` among its free symbols, we leave it alone.
    relations : list of sympy.Relational
        A list of known inequalities among your constant symbols, e.g. [L > a, b > a].
        Supported types are Gt (>) / Ge (≥) / Lt (<) / Le (≤).  The code uses these
        to decide the sign of `(f - knot)` whenever `f` is truly a constant (no `x`).

    Returns
    -------
    sympy.Basic
        A new expression in which all “collapsible” SingularityFunction nodes have
        been replaced by either `0`, `1`, or `(f - knot)**n` as appropriate.
    """
    def _collapse_sf(node):
        # Only process SingularityFunction nodes
        if not isinstance(node, SingularityFunction):
            return node

        f, knot, n = node.args

        # 1) If the “variable” appears in f, leave this singularity alone
        if var in f.free_symbols:
            return node

        # 2) Otherwise, compute diff = f - knot and simplify
        diff = simplify(f - knot)

        # 3) If diff.is_positive or diff.is_negative is already decided by Sympy, use that:
        if diff.is_positive:
            # f > knot
            return (1 if n == 0 else diff**n)
        if diff.is_negative:
            # f < knot
            return 0  # for n >= 0

        # 4) If still undecided, see if any of our user‐supplied relations matches “diff”
        for rel in relations:
            if not isinstance(rel, Relational):
                continue
            u, v = rel.lhs, rel.rhs
            target = simplify(u - v)

            # Case A: diff ==  (u - v)
            if simplify(diff - target) == 0:
                # Now dispatch based on whether rel is >, ≥, <, or ≤
                if isinstance(rel, Gt):   # u > v  => (u - v) > 0
                    return (1 if n == 0 else diff**n)
                if isinstance(rel, Ge):   # u ≥ v => (u - v) ≥ 0
                    if n == 0:
                        return 1
                    else:
                        return diff**n
                if isinstance(rel, Lt):   # u < v  => (u - v) < 0
                    return 0
                if isinstance(rel, Le):   # u ≤ v => (u - v) ≤ 0
                    if diff == 0 and n == 0:
                        return 1
                    return 0

            # Case B: diff == -(u - v)  (i.e. diff = v - u)
            if simplify(diff + target) == 0:
                if isinstance(rel, Gt):   # u > v  => (v - u) < 0
                    return 0
                if isinstance(rel, Ge):   # u ≥ v => (v - u) ≤ 0
                    if diff == 0 and n == 0:
                        return 1
                    return 0
                if isinstance(rel, Lt):   # u < v  => (v - u) > 0
                    return (1 if n == 0 else diff**n)
                if isinstance(rel, Le):   # u ≤ v => (v - u) ≥ 0
                    if n == 0:
                        return 1
                    return diff**n

        # 5) If we still can’t decide sign, leave it alone
        return node

    # Traverse the entire expression and replace every SingularityFunction by _collapse_sf(...)
    return expr.replace(lambda e: isinstance(e, SingularityFunction), _collapse_sf)

def format_label(text: str) -> str:
    """
    Convert a plain‐text label into Plotly‐HTML with:
      1) simple fractions     a/b  →  <sup>a</sup>/<sub>b</sub>
      2) underscores          x_y  →  x<sub>y</sub>
      3) named Greek letters  phi, theta, etc.  → their Unicode (φ, θ, …)
    
    NOTE on ordering:
    - We do the fraction‐replacement first.
    - Then we do the underscore‐to‐<sub>…</sub> replacement.
    - Finally, with all HTML tags in place, we do a simple string‐replace
      for any Greek names (like “phi”) → their Unicode “φ”. This way,
      if the subscript token is “pphi”, it becomes “p<sub>pphi</sub>” first,
      and only then we replace “phi” → “ϕ” inside that subscript.
    """
    text = text.replace('{','').replace('}','')

    # 1) Fractions: replace a/b ⇒ <sup>a</sup>/<sub>b</sub>
    def _frac(match):
        num, den = match.group(1), match.group(2)
        return f"<sup>{num}</sup>/<sub>{den}</sub>"

    text = re.sub(
        r'\b([A-Za-z0-9]+)\/([A-Za-z0-9]+)\b',
        _frac,
        text
    )

    # 2) Subscripts: replace every occurrence of "_token" ⇒ "<sub>token</sub>"
    text = re.sub(r'_(\w+)', r'<sub>\1</sub>', text)

    # 3) Named‐Greek → Unicode
    #    (only lowercase names are handled here; add uppercase if needed)
    greek_map = {
        "alpha":   "α",
        "beta":    "β",
        "gamma":   "γ",
        "delta":   "δ",
        "epsilon": "ε",
        "zeta":    "ζ",
        "theta":   "θ",
        "iota":    "ι",
        "kappa":   "κ",
        "lambda":  "λ",
        "mu":      "μ",
        "nu":      "ν",
        "xi":      "ξ",
        "omicron": "ο",
        "pi":      "π",
        "rho":     "ρ",
        "sigma":   "σ",
        "tau":     "τ",
        "upsilon": "υ",
        "phi":     "ϕ",
        "chi":     "χ",
        "psi":     "ψ",
        "omega":   "ω",
    }

    # Simply replace ANY occurrence of “phi” (etc.) with its Unicode.
    # Because we already wrapped things in <sub>…</sub> or <sup>…</sup>, 
    # this will correctly transform “pphi” → “pϕ” inside that tag.
    for name, uni in greek_map.items():
        text = text.replace(name, uni)

    return text

def format_subs(label: str) -> str:
    label = label.replace('{','').replace('}','')
    a = label.split('_')
    if len(a) > 1:
        total_subs = len(a) - 1
        label = '<sub>'.join(a)
        for i in range(total_subs):
            label += '</sub>'
    return label

def latex_with_threshold(expr, **printer_kwargs):
    """
    Render expr to LaTeX, using scientific notation only if
    exponent > sci_max or < sci_min (defaults +5 / −2).
    Works on sympy expressions *or* plain Python floats/ints.
    """
    # 1) Turn Python native types into SymPy objects
    expr = sympify(expr)    
    # 2) Create and use our custom printer
    printer = ThresholdLatexPrinter(**printer_kwargs)
    return printer.doprint(expr)

def fig_to_rotated_img(fig: go.Figure, rotation_angle: int = 0) -> str:
    """
    Convert a Plotly figure to a rotated PNG image (in-memory) and return it as an <img> HTML element.

    Parameters:
    - fig: Plotly Figure object
    - rotation_angle: degrees to rotate clockwise (positive) or counterclockwise (negative)

    Returns:
    - HTML string for a Data URI <img> tag with the rotated image
    """
    # Export Plotly figure to PNG bytes
    png_bytes = fig.to_image(format="png")

    # Rotate image in memory using Pillow
    with Image.open(io.BytesIO(png_bytes)) as im:
        # Pillow uses counterclockwise for positive angle; invert if needed
        rotated = im.rotate(-rotation_angle, expand=True)
        buffer = io.BytesIO()
        rotated.save(buffer, format="PNG")
        buffer.seek(0)
        b64_str = base64.b64encode(buffer.read()).decode("utf-8")

    # Form Data URI and wrap in an <img> tag
    data_uri = f"data:image/png;base64,{b64_str}"
    img_html = f'<img id="model" src="{data_uri}" alt="Rotated Plot" />'
    return img_html

def parse_and_update_symbols(expr_str: str, existing: list[str]) -> list[str]:
    """
    Parse expr_str with sympify, re‑using any names in `existing`,
    and return a list of any new free symbols (as strings).
    """
    # 1) tell sympify about your existing names
    locals_dict = { name: Symbol(name) for name in existing }
    
    # 2) parse — unknown names become fresh Symbols
    expr = sympify(expr_str, locals=locals_dict)
    
    # 3) collect all symbols seen in the result
    all_names = { str(s) for s in expr.free_symbols }
    
    # 4) subtract the ones you already had
    new_names = all_names - set(existing)
    
    return list(new_names)

def is_only_one_L_and_numbers(expr_str: str, primary: str = "L") -> bool:
    """
    Return True if `expr_str` parses to a SymPy expression whose only free symbol
    is `primary` (e.g. "L") and that symbol appears exactly once in the entire tree.
    Otherwise return False.

    Examples:
        is_only_one_L_and_numbers("L - 5")     -> True
        is_only_one_L_and_numbers("L - L/2")   -> False   (L appears twice)
        is_only_one_L_and_numbers("L - a")     -> False   (a is an extra symbol)
        is_only_one_L_and_numbers("2*L + 3")   -> False   (L appears once? Actually count is 1—True)
        is_only_one_L_and_numbers("2*L")       -> True
        is_only_one_L_and_numbers("5")         -> False   (no L at all)
    """
    # 1) Create a SymPy Symbol for “L” (or whatever `primary` is)
    L = symbols(primary)

    # 2) Parse the string into a SymPy expression, telling sympify that
    #    any occurrence of the name `primary` is the symbol L.  If the string
    #    contains any other letter‐strings (like “a” or “x”), those become new Symbols
    #    and end up in expr.free_symbols.
    try:
        expr = sympify(expr_str, locals={primary: L})
    except Exception:
        # If parsing fails entirely, we treat that as “not valid”
        return False

    # 3) Check which symbols actually appear:
    free_syms = expr.free_symbols
    #    We want exactly {L} (no extras).  If free_syms is empty or contains anything
    #    other than L, return False.
    if free_syms != {L}:
        return False

    # 4) Count how many times L appears in the expression tree.
    #    expr.count(L) returns the number of sub‐expressions equal to L.
    #    For example:
    #       - sympify("L - 5").count(L)  == 1
    #       - sympify("L - L/2").count(L) == 2
    #       - sympify("2*L").count(L) == 1
    cnt = expr.count(L)
    return (cnt == 1)

def system_to_latex(*equations):
    if not equations:
        return ""
    # Join the LaTeX representations of each equation with line breaks
    equations_latex = r" \\ ".join(latex_with_threshold(eq) for eq in equations)
    # Wrap the equations in a LaTeX system environment
    return r"\left\{\begin{matrix}" + equations_latex + r"\end{matrix}\right."

def arrow_inline_no_braces(eq1, eq2, arrow="\\rightarrow"):
    """
    Gera código LaTeX inline: eq1 arrow eq2, sem chaves.
    - eq1, eq2: objetos sympy Eq()
    - arrow: símbolo da seta, ex: "\\rightarrow", "\\Rightarrow"
    """
    l1 = latex_with_threshold(eq1)
    l2 = latex_with_threshold(eq2)
    return rf"{l1} \; {arrow} \; {l2}"
