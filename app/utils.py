from sympy import Float, sympify, Symbol, Piecewise, Eq, simplify
from sympy.assumptions.ask import ask
from sympy.assumptions import Q
from sympy.functions.special.singularity_functions import SingularityFunction
from sympy.printing.latex import LatexPrinter
from typing import List, Union
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
        # 1) Turn the SymPy Float into a native Python float
        pyf = float(f)

        # 2) Use Python’s "{:e}" formatting to get "mantissa e±exp"
        sci_str = "{:.15e}".format(pyf)               # e.g. "1.234567890123450e+06" :contentReference[oaicite:0]{index=0}
        mant_str, exp_str = sci_str.split("e")
        e = int(exp_str)

        # 3) Strip trailing zeros & “.” from the mantissa
        mant_str = mant_str.rstrip("0").rstrip(".")

        # 4) If exponent outside [sci_min, sci_max], emit in m×10^e form
        if e > self.sci_max or e < self.sci_min:
            return rf"{mant_str}\times10^{{{e}}}"

        # 5) Otherwise defer to SymPy’s default float printer
        return super()._print_Float(f)
    
    def _print_SingularityFunction(self, expr):
        xx, aa, nn = expr.args
        return r'{\left\langle %s - %s \right\rangle}^{%s}' % (
            self._print(xx),  # x
            self._print(aa),  # a
            self._print(nn),  # n
        )

Item = List[Union[str, float, int]]

def add_item(collection: List[Item], new_item: Item) -> bool:
    """
    Add `new_item` to `collection` if there is no existing item in `collection`
    whose first three elements match new_item's first three elements.
    
    Returns True if new_item was added, False otherwise.
    """
    # Compare only the first three elements
    key_new = tuple(new_item[:3])
    for item in collection:
        if tuple(item[:3]) == key_new:
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

def eval_all_singularities(expr_eq: Eq, sub_symbol, sub_value):
    """
    For Eq(lhs, rhs) where rhs may be a sum of SingularityFunction(...) + const,
    substitute sub_symbol → sub_value *inside* each SF, force its Piecewise
    evaluation, and simplify the final result.
    """
    def _eval_one(sf: SingularityFunction):
        # 1) rewrite this single SF → Piecewise, then substitute
        pw = sf.rewrite(Piecewise).subs(sub_symbol, sub_value)
        # 2) grab the branch‐condition (assumed on index 0)
        cond = pw.args[0][1]
        # 3) decide which branch holds under sub_value > 0
        if ask(cond, Q.is_true(sub_value > 0)):
            return pw.args[0][0]
        else:
            return pw.args[1][0]

    # 4) apply to every SingularityFunction in the RHS
    new_rhs = expr_eq.rhs.replace(
        lambda e: isinstance(e, SingularityFunction),
        lambda e: _eval_one(e)
    )
    # 5) simplify any remaining arithmetic
    new_rhs = simplify(new_rhs)
    # 6) rebuild the equality
    return Eq(expr_eq.lhs, new_rhs)

def format_label(text: str) -> str:
    """
    Convert a plain‐text label into Plotly‐HTML with:
      * simple fractions a/b → <sup>a</sup>/<sub>b</sub>
      * underscores  x_y → x<sub>y</sub> (only the token after each `_`)
    """
    # 1) Fractions: <sup>num</sup>/<sub>den</sub>
    def _frac(m):
        num, den = m.group(1), m.group(2)
        return f"<sup>{num}</sup>/<sub>{den}</sub>"

    text = re.sub(
        r'\b([A-Za-z0-9]+)\/([A-Za-z0-9]+)\b',
        _frac,
        text
    )

    # 2) Subscripts: replace every _token with <sub>token</sub>
    #    \w+ matches letters, digits or underscore — adjust if you need hyphens etc.
    text = re.sub(r'_(\w+)', r'<sub>\1</sub>', text)

    return text

def format_subs(label: str) -> str:
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
    img_html = f'<img src="{data_uri}" alt="Rotated Plot" />'
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