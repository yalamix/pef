import plotly.graph_objects as go
import numpy as np

def plot_rectangle(width, height,
                                          rect_line_color='black', 
                                          rect_fillcolor='rgba(0,0,0,0)', 
                                          rect_line_width=3,
                                          arrow_length=1,
                                          arrow_line_width=1):
    """
    Draws a rectangle starting at (0,0) with given width and height, and
    adds fixed-length x and y axis arrows (both of length arrow_length, default 1)
    that start at (0,0). The arrow lines are thinner, drawn with semi-transparent
    black, and the "x" and "y" labels are placed near the arrow tips with shifts
    to avoid overlapping the arrow lines.
    
    Parameters:
      width (float): Width of the rectangle.
      height (float): Height of the rectangle.
      rect_line_color (str): Rectangle edge color.
      rect_fillcolor (str): Rectangle fill color.
      rect_line_width (float): Thickness of rectangle edge.
      arrow_length (float): Length of both x and y axis arrows (fixed to 1 unit by default).
      arrow_line_width (float): Thickness of the arrow lines.
      
    Returns:
      fig (go.Figure): The Plotly figure.
    """
    # Rectangle: from (0,0) to (width, height)
    x1 = width
    y1 = height

    fig = go.Figure()

    # Draw the rectangle
    fig.add_shape(
        type="rect",
        x0=0, y0=0, x1=x1, y1=y1,
        line=dict(color=rect_line_color, width=rect_line_width),
        fillcolor=rect_fillcolor
    )

    # Update axes to include both the rectangle and the fixed-length arrows.
    x_min = -width * 0.2
    x_max = max(width, arrow_length) + width * 0.2
    y_min = -height * 0.2
    y_max = max(height, arrow_length) + height * 0.2
    fig.update_xaxes(range=[x_min, x_max])
    fig.update_yaxes(range=[y_min, y_max])
    
    # Enforce equal scaling so that one unit in x equals one unit in y.
    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    
    fig.update_layout(
        # title="Barra",
        template="plotly_white",
        width=960,
        height=720
    )
    
    return fig

def add_vector(fig, tail, head, label=None, arrow_color='black', arrow_width=2, arrow_head=3, xanchor="auto", yanchor="top"):
    """
    Add a vector as an arrow with a text label to a Plotly figure.
    
    Parameters:
      fig (go.Figure): Plotly figure to which the arrow is added.
      tail (tuple): (x, y) coordinates where the arrow tail (start) is.
      head (tuple): (x, y) coordinates where the arrow head (end) is.
      label (str): Text to display near the arrow.
      arrow_color (str): Color of the arrow. Default is 'black'.
      arrow_width (int): Width of the arrow. Default is 2.
      arrow_head (int): Style number of the arrow head. Default is 3.
      
    Returns:
      go.Figure: Updated Plotly figure with the vector added.
    """
    # Add an annotation that creates an arrow from tail to head
    fig.add_annotation(
        x=head[0],  # arrow head x-coordinate
        y=head[1],  # arrow head y-coordinate
        ax=tail[0],  # arrow tail x-coordinate
        ay=tail[1],  # arrow tail y-coordinate
        text=label,  # text label to display near the arrow
        showarrow=True,
        arrowhead=arrow_head,
        arrowwidth=arrow_width,
        arrowcolor=arrow_color,
        axref="x",
        ayref="y",
        xref="x",
        yref="y",
        xanchor=xanchor,
        yanchor=yanchor,
        font=dict(color=arrow_color, size=16)
    )
    return fig

def add_axis_arrows(fig):
    """
    Adds x and y axis arrows.
    
    Parameters:
      fig         : plotly.graph_objects.Figure object to add traces to.
    
    Returns:
      The updated Plotly figure.
    """    
    fig = add_vector(fig, (0, 0), (0.25, 0), 'x', "rgba(0,0,0,0.5)")
    fig = add_vector(fig, (0, 0), (0, 0.25), 'y', "rgba(0,0,0,0.5)", xanchor="right", yanchor="bottom")
    return fig

def add_curve_with_y_cutoff_fill(
    fig, x, y, y_cutoff,
    annotation_text=None,
    annotation_offset=None,
    curve_color='rgba(0,100,80,1)', 
    fill_color='rgba(0,100,80,0.2)', 
    tol=1e-8,
    up=True,
    side=False,
    double=False
):
    """
    Adds a curve to the figure and fills the area between the curve and a horizontal
    line at y = y_cutoff for regions where the curve is above the cutoff.
    Optionally, an annotation is placed above the curve.
    
    Parameters:
      fig               : Plotly Figure object.
      x                 : List or array of x values.
      y                 : List or array of y values.
      y_cutoff          : The y value at which to cut off the fill.
      annotation_text   : (Optional) Text to annotate above the curve.
      annotation_offset : (Optional) Vertical offset to place the annotation above the curve.
                          If not provided, it is computed as 10% of the y-range.
      curve_color       : Color for the curve line.
      fill_color        : Color for the filled area under the curve.
      tol               : Tolerance for detecting degenerate segments.
      up               : Direction of arrows.
      
    Returns:
      Updated Plotly Figure.
    """
    
    # Add the full curve trace with legend hidden.
    fig.add_trace(go.Scatter(
        x=x,
        y=y,
        mode='lines',
        showlegend=False,
        line=dict(color=curve_color)
    ))
    
    # Helper: linear interpolation for crossing the y_cutoff.
    def interp_cross(x0, y0, x1, y1, y_cut):
        t = (y_cut - y0) / (y1 - y0)
        return x0 + t * (x1 - x0)
    
    segments = []  # To store contiguous segments where the curve is above y_cutoff.
    current_segment = []
    N = len(x)
    
    for i in range(N-1):
        p0 = (x[i], y[i])
        p1 = (x[i+1], y[i+1])
        above0 = p0[1] >= y_cutoff
        above1 = p1[1] >= y_cutoff
        
        # If p0 is above the cutoff, add it.
        if above0:
            if not current_segment or current_segment[-1] != p0:
                current_segment.append(p0)
        
        # Check if the segment between p0 and p1 crosses y_cutoff.
        if above0 != above1:
            xi = interp_cross(p0[0], p0[1], p1[0], p1[1], y_cutoff)
            crossing = (xi, y_cutoff)
            
            if above0:
                # Curve goes from above to below: add crossing and finish segment.
                current_segment.append(crossing)
                segments.append(current_segment)
                current_segment = []
            else:
                # Curve goes from below to above: start new segment with crossing.
                current_segment = [crossing]
    
    # If an active segment remains at the end, add the last point if above the cutoff.
    if current_segment:
        if y[-1] >= y_cutoff:
            current_segment.append((x[-1], y[-1]))
        segments.append(current_segment)
    
    # Add filled areas for each non-degenerate segment.
    for seg in segments:
        if len(seg) < 2:
            continue
        
        # Skip degenerate segments (e.g., constant at the cutoff).
        if all(np.abs(pt[1] - y_cutoff) < tol for pt in seg):
            continue
        
        seg_x = [pt[0] for pt in seg]
        seg_y = [pt[1] for pt in seg]
        
        # Create polygon by following the curve and closing along y = y_cutoff.
        poly_x = seg_x + [seg_x[-1]]
        poly_y = seg_y + [y_cutoff]
        poly_x = [seg_x[0]] + poly_x
        poly_y = [y_cutoff] + poly_y
        
        fig.add_trace(go.Scatter(
            x=poly_x,
            y=poly_y,
            mode='lines',
            showlegend=False,
            fill='toself',
            fillcolor=fill_color,
            line=dict(color='rgba(0,0,0,0)')  # Hide boundary
        ))
    
    # Add annotation above the curve if annotation_text is provided.
    if annotation_text is not None:
        # Determine a reasonable position:
        annotation_x = np.mean(x)
        max_y = np.max(y)
        min_y = np.min(y)
        # If the curve is constant, choose a small fixed offset.
        if np.isclose(max_y, min_y):
            default_offset = 0.1 * (abs(max_y) if max_y != 0 else 1)
        else:
            default_offset = 0.1 * (max_y - min_y)
        offset = annotation_offset if annotation_offset is not None else default_offset
        annotation_y = max_y + offset
        
        fig.add_annotation(
            x=annotation_x,
            y=annotation_y,
            text=annotation_text,
            showarrow=False,
            font=dict(size=16, color=curve_color),
            align="center"
        )
    
    # Add arrows
    a = int((np.max(x) - np.min(x))/0.2)
    b_size = len(x)//a    

    if side:
        # Normal force and twisting moment
        a = int((2 if double else 1) *(np.max(x) - np.min(x))/0.2)
        b_size = len(x)//a    
        x_size = (x[b_size] - x[0])/2
        for i in range(a):
            if up:
                if (y[i * b_size] + y[(i+1) * b_size if (i+1) * b_size < len(x) else i * b_size])/2 > y_cutoff + 0.1 and (i % 2 == 0 if double else True):
                    fig = add_vector(fig, (x[i * b_size], y_cutoff + 0.1), (x[(i + 1) * b_size], y_cutoff + 0.1), arrow_color=curve_color, arrow_width=1.5)
                    if double:
                        fig = add_vector(fig, (x[i * b_size] + x_size, y_cutoff + 0.1), (x[(i + 1) * b_size] + x_size, y_cutoff + 0.1), arrow_color=curve_color, arrow_width=1.5)
            else:
                if (y[i * b_size] + y[(i+1) * b_size if (i+1) * b_size < len(x) else i * b_size])/2 > y_cutoff + 0.1 and (i % 2 == 0 if double else True):
                    fig = add_vector(fig, (x[(i + 1) * b_size], y_cutoff + 0.1), (x[i * b_size], y_cutoff + 0.1), arrow_color=curve_color, arrow_width=1.5)            
                    if double:
                        fig = add_vector(fig, (x[(i + 1) * b_size] + x_size, y_cutoff + 0.1), (x[i * b_size] + x_size, y_cutoff + 0.1), arrow_color=curve_color, arrow_width=1.5)            
    else:
        for i in range(a):
            if up:
                fig = add_vector(fig, (x[i * b_size], y_cutoff), (x[i * b_size], y[i * b_size]), arrow_color=curve_color)
            else:
                fig = add_vector(fig, (x[i * b_size], y[i * b_size]), (x[i * b_size], y_cutoff), arrow_color=curve_color)
        if up:
            fig = add_vector(fig, (x[-1], y_cutoff), (x[-1], y[-1]), arrow_color=curve_color)
        else:
            fig = add_vector(fig, (x[-1], y[-1]), (x[-1], y_cutoff), arrow_color=curve_color)                
    
    return fig

def add_cantilever(fig, rec_height, x = 0):
    """
    Adds a cantilever to the beam.
    
    Parameters:
      fig         : plotly.graph_objects.Figure object to add traces to.
      rec_height          : The height of the beam.
      x          : x position of the cantilever.
    
    Returns:
      The updated Plotly figure.
    """    
    c_size = rec_height * 3
    f_size = c_size/10

    fig.add_shape(type="line",
        xref="x", yref="y",
        x0=x, y0=-c_size/2 + rec_height/2, x1=x, y1=+c_size/2 + rec_height/2,
        line=dict(
            color="black",
            width=3,
        ),
    )

    for i in range(9):
        fig.add_shape(type="line",
            xref="x", yref="y",
            x0=x - f_size if x == 0 else x + f_size, y0=+c_size/2 + rec_height/2   - i * f_size, x1=x, y1=+c_size/2 + rec_height/2  - (i+1) * f_size,
            line=dict(
                color="black",
                width=3,
            ),
        )

    return fig

def add_hinge(fig, rec_height, x = 0):
    """
    Adds a hinge to the beam.
    
    Parameters:
      fig         : plotly.graph_objects.Figure object to add traces to.
      rec_height          : The height of the beam.
      x          : x position of the hinge.
    
    Returns:
      The updated Plotly figure.
    """    
    r = rec_height * 0.95
    fig.add_shape(type="circle",
        xref="x", yref="y",
        x0=x - r/2, y0=rec_height/2 - r/2, x1=x + r/2, y1=rec_height/2 + r/2,
        line_color="DarkOrange",
    )

    return fig

def add_fixed_support(fig, rec_height, x = 0):
    """
    Adds a fixed hinged support to the beam.
    
    Parameters:
      fig         : plotly.graph_objects.Figure object to add traces to.
      rec_height          : The height of the beam.
      x          : x position of the support.
    
    Returns:
      The updated Plotly figure.
    """    
    c_size = rec_height * 2
    f_size = c_size/8
    fig.add_shape(type="line",
        xref="x", yref="y",
        x0=x - rec_height/2, y0=-rec_height, x1=x, y1=0,
        line=dict(
            color="black",
            width=3,
        ),
    )
    fig.add_shape(type="line",
        xref="x", yref="y",
        x0=x + rec_height/2, y0=-rec_height, x1=x, y1=0,
        line=dict(
            color="black",
            width=3,
        ),
    )
    fig.add_shape(type="line",
        xref="x", yref="y",
        x0=x - rec_height, y0=-rec_height, x1=x + rec_height, y1=-rec_height,
        line=dict(
            color="black",
            width=3,
        ),
    )

    for i in range(1,7):
        fig.add_shape(type="line",
            xref="x", yref="y",
            x0=x - rec_height + i * f_size, y0=-rec_height, x1=x - rec_height + (i + 1) * f_size, y1=-rec_height - f_size,
            line=dict(
                color="black",
                width=3,
            ),
        )

    return fig

def add_mobile_support(fig, rec_height, x = 0):
    """
    Adds a mobile hinged support to the beam.
    
    Parameters:
      fig         : plotly.graph_objects.Figure object to add traces to.
      rec_height          : The height of the beam.
      x          : x position of the support.
    
    Returns:
      The updated Plotly figure.
    """    
    c_size = rec_height * 2
    f_size = c_size/8
    fig.add_shape(type="line",
        xref="x", yref="y",
        x0=x - rec_height/2, y0=-rec_height + f_size, x1=x, y1=0,
        line=dict(
            color="black",
            width=3,
        ),
    )
    fig.add_shape(type="line",
        xref="x", yref="y",
        x0=x + rec_height/2, y0=-rec_height + f_size, x1=x, y1=0,
        line=dict(
            color="black",
            width=3,
        ),
    )
    fig.add_shape(type="line",
        xref="x", yref="y",
        x0=x - rec_height/2, y0=-rec_height + f_size, x1=x + rec_height/2, y1=-rec_height + f_size,
        line=dict(
            color="black",
            width=3,
        ),
    )    
    fig.add_shape(type="circle",
        xref="x", yref="y",
        x0=x - rec_height/2, y0=-rec_height, x1=x - rec_height/2 + f_size, y1=-rec_height + f_size,
        line_color="black",
    )
    fig.add_shape(type="circle",
        xref="x", yref="y",
        x0=x + rec_height/2 -f_size, y0=-rec_height, x1=x + rec_height/2, y1=-rec_height + f_size,
        line_color="black",
    )    
    fig.add_shape(type="line",
        xref="x", yref="y",
        x0=x - rec_height, y0=-rec_height, x1=x + rec_height, y1=-rec_height,
        line=dict(
            color="black",
            width=3,
        ),
    )

    for i in range(1,7):
        fig.add_shape(type="line",
            xref="x", yref="y",
            x0=x - rec_height + i * f_size, y0=-rec_height, x1=x - rec_height + (i + 1) * f_size, y1=-rec_height - f_size,
            line=dict(
                color="black",
                width=3,
            ),
        )

    return fig

def add_roller_support(fig, rec_height, x=0, roller_radius_frac=0.08, bracket_width_frac=0.2, bracket_overhang_frac=0.05):
    """
    Adds a roller support (two small circles + surrounding bracket) centered at x.
    
    Parameters:
      fig                   : plotly.graph_objects.Figure to add shapes to.
      rec_height            : the full height of your beam (y-extent).
      x                     : the beam-wise center of the support.
      roller_radius_frac    : radius of each roller as a fraction of rec_height.
      bracket_width_frac    : half-width of the bracket arms, as fraction of rec_height.
      bracket_overhang_frac : extra bracket length above top roller & below bottom roller,
                              as fraction of rec_height.
    
    Returns:
      The updated Figure.
    """
    # absolute sizes
    r = rec_height * roller_radius_frac
    half_bw = rec_height * bracket_width_frac
    overhang = rec_height * bracket_overhang_frac

    # vertical positions of the two rollers (stacked)
    y_center = rec_height / 2
    y1 = y_center + r * 2
    y2 = y_center - r * 2

    # bracket top & bottom
    y_top = y1 + r + overhang - r/2
    y_bot = y2 - r - overhang + r/2

    # x-coords of left & right bracket arms
    xL = x - half_bw/2
    xR = x + half_bw/2

    # draw the two rollers
    for yc in (y1, y2):
        fig.add_shape(
            type="circle",
            xref="x", yref="y",
            x0=x - r, y0=yc - r,
            x1=x + r, y1=yc + r,
            line_color="black",
            fillcolor="lightgray",
        )

    # draw left vertical bracket arm
    fig.add_shape(
        type="line",
        xref="x", yref="y",
        x0=xL, y0=y_bot - y_center,
        x1=xL, y1=y_top + y_center,
        line_color="black",
        line_width=2,
    )
    # little horizontal ticks on left arm
    fig.add_shape(type="line", xref="x", yref="y",
                  x0=xL,             y0=y_top + y_center,
                  x1=xL + 2 * half_bw, y1=y_top + y_center,
                  line_color="black", line_width=2)
    fig.add_shape(type="line", xref="x", yref="y",
                  x0=xL,             y0=y_bot - y_center,
                  x1=xL + 2 * half_bw, y1=y_bot - y_center,
                  line_color="black", line_width=2)

    # little vertical ticks on left arm
    fig.add_shape(
        type="line",
        xref="x", yref="y",
        x0=xL + 2 * half_bw, y0=y_top + y_center - half_bw,
        x1=xL + 2 * half_bw, y1=y_top + y_center,
        line_color="black",
        line_width=2,
    )
    fig.add_shape(
        type="line",
        xref="x", yref="y",
        x0=xL + 2 * half_bw, y0=y_bot - y_center,
        x1=xL + 2 * half_bw, y1=y_bot - y_center + half_bw,
        line_color="black",
        line_width=2,
    )    

    # draw right vertical bracket arm
    fig.add_shape(
        type="line",
        xref="x", yref="y",
        x0=xR, y0=y_bot - half_bw,
        x1=xR, y1=y_top + half_bw,
        line_color="black",
        line_width=2,
    )

    return fig

def add_semicircle_arrow(
    fig,
    x_center,
    y_center,
    radius=0.5,
    orientation='counterclockwise',
    arrow_size=0.1,
    arc_points=100,
    color='black',
    start_angle=np.radians(270),
    label=None,
    label_offset=0.3,
    label_font=dict(color="black", size=14)
):
    """
    Adds a semicircular arrow to the given Plotly figure, with an optional label.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The figure to which the semicircular arrow will be added.
    x_center, y_center : float
        Center of the semicircle.
    radius : float
        Radius of the semicircle.
    orientation : str
        'counterclockwise' or 'clockwise'.
    arrow_size : float
        Length of the arrowhead lines.
    arc_points : int
        Number of points to use for drawing the arc (smoothness).
    color : str
        Color of the semicircular arrow.
    start_angle : float
        Starting angle (in radians) of the semicircle, measured from the positive x-axis.
    label : str, optional
        If provided, adds this text label above the semicircle.
    label_offset : float, optional
        Radial offset added to the arc's radius for label placement (default is 0.3).
    label_font : dict, optional
        Font styling for the label.
    """
    # Determine angle range for a semicircle
    if orientation == 'counterclockwise':
        angles = np.linspace(start_angle, start_angle + np.pi, arc_points)
    else:
        angles = np.linspace(start_angle, start_angle - np.pi, arc_points)

    # Parametric semicircle
    x_arc = x_center + radius * np.cos(angles)
    y_arc = y_center + radius * np.sin(angles)

    # Add the arc as a line
    fig.add_trace(
        go.Scatter(
            x=x_arc,
            y=y_arc,
            mode='lines',
            line=dict(color=color),
            showlegend=False
        )
    )

    # Compute the tangent direction at the end of the arc
    theta_end = angles[-1]
    dx = -radius * np.sin(theta_end)
    dy =  radius * np.cos(theta_end)
    # Normalize the tangent vector
    norm = np.sqrt(dx**2 + dy**2)
    dx, dy = dx/norm, dy/norm

    # The arrow tip is the last point on the arc
    x_tip = x_arc[-1]
    y_tip = y_arc[-1]

    # Create a small “V” shape for the arrowhead.
    arrow_half_angle = np.radians(20)  # ~20° spread

    def rotate(vx, vy, alpha):
        """Rotate vector (vx, vy) by angle alpha (radians)."""
        rx = vx * np.cos(alpha) - vy * np.sin(alpha)
        ry = vx * np.sin(alpha) + vy * np.cos(alpha)
        return rx, ry

    # Get the two directions for the arrowhead
    dx_left, dy_left = rotate(dx, dy,  arrow_half_angle)
    dx_right, dy_right = rotate(dx, dy, -arrow_half_angle)

    # Calculate endpoints for the arrowhead lines
    if orientation == 'counterclockwise':
        x_left  = x_tip - arrow_size * dx_left
        y_left  = y_tip - arrow_size * dy_left
        x_right = x_tip - arrow_size * dx_right
        y_right = y_tip - arrow_size * dy_right
    else:
        x_left  = x_tip + arrow_size * dx_left
        y_left  = y_tip + arrow_size * dy_left
        x_right = x_tip + arrow_size * dx_right
        y_right = y_tip + arrow_size * dy_right        

    # Add arrowhead traces
    fig.add_trace(
        go.Scatter(
            x=[x_tip, x_left],
            y=[y_tip, y_left],
            mode='lines',
            line=dict(color=color),
            showlegend=False
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[x_tip, x_right],
            y=[y_tip, y_right],
            mode='lines',
            line=dict(color=color),
            showlegend=False
        )
    )

    # If label is provided, add it as an annotation above the arc.
    if label is not None:
        # Compute the midpoint angle for the arc.
        if orientation == 'counterclockwise':
            mid_angle = start_angle + np.pi/2
        else:
            mid_angle = start_angle - np.pi/2

        # Position the label at the arc's midpoint, then push it outward by label_offset.
        x_label = x_center + (radius + label_offset) * np.cos(mid_angle)
        y_label = y_center + (radius + label_offset) * np.sin(mid_angle)
        fig.add_annotation(
            x=x_label,
            y=y_label,
            text=label,
            showarrow=False,
            font=label_font,
            xanchor="center"
        )

    return fig
