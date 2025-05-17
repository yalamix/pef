from sympy import Float, Rational, sympify
from sympy.printing.latex import LatexPrinter
from PIL import Image
import plotly.graph_objects as go
import io
import base64

class ThresholdLatexPrinter(LatexPrinter):
    def __init__(self, *args, sci_min=-2, sci_max=5, **kwargs):
        super().__init__(*args, **kwargs)
        self.sci_min = sci_min
        self.sci_max = sci_max

    def _print_Float(self, f: Float) -> str:
        # 1) Recover exact numerator/denominator of the Python float
        py_num, py_den = float(f).as_integer_ratio()
        rat = Rational(py_num, py_den)

        # 2) Compute mantissa (m) and exponent (e) in base-10:
        if rat.q == 1:
            # rat is an integer; factor out powers of 10 manually
            m_int = rat.p
            e = 0
            while m_int % 10 == 0 and m_int != 0:
                m_int //= 10
                e += 1
            m = m_int
        else:
            # a true rational: use as_mant_exp()
            m, e = rat.as_mant_exp()

        # 3) Decide formatting based on thresholds
        if e > self.sci_max or e < self.sci_min:
            # format mantissa, strip “.0” if any
            m_str = str(m)
            if m_str.endswith('.0'):
                m_str = m_str[:-2]
            return rf"{m_str}\times10^{{{e}}}"
        # 4) Otherwise defer to SymPy’s default float‐printer
        return super()._print_Float(f)

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
