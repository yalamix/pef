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
