import matplotlib as mpl

# --- IEEE Figure Dimensions (in inches) ---
IEEE_SINGLE_COLUMN_WIDTH = 3.5  # inches
IEEE_DOUBLE_COLUMN_WIDTH = 7.16  # inches

# --- IEEE Font Sizes (in points) ---
# Main body text in IEEE is typically 10pt.
# Figure text should be comparable or slightly smaller.
FONT_SIZE_SMALL = 8  # For tick labels, legend entries
FONT_SIZE_MEDIUM = 9  # For axis labels, inline text
FONT_SIZE_LARGE = 10  # For sub-plot titles


def apply_ieee_mpl_settings(column_width="single"):
    if column_width == "single":
        fig_width = IEEE_SINGLE_COLUMN_WIDTH
    elif column_width == "double":
        fig_width = IEEE_DOUBLE_COLUMN_WIDTH
    else:
        raise ValueError("column_width must be 'single' or 'double'")

    fig_height = fig_width / 1.618

    mpl.rcParams["mathtext.fontset"] = "cm"  # Computer Modern for math text
    mpl.rcParams["font.family"] = "STIXGeneral"

    # --- Figure Dimensions ---
    mpl.rcParams["figure.figsize"] = (fig_width, fig_height)
    mpl.rcParams["figure.dpi"] = (
        300  # 300 for standard print resolution (or 600 for line art)
    )

    # --- Font Sizes ---
    mpl.rcParams["font.size"] = FONT_SIZE_MEDIUM
    mpl.rcParams["axes.labelsize"] = FONT_SIZE_SMALL
    mpl.rcParams["axes.titlesize"] = FONT_SIZE_MEDIUM
    mpl.rcParams["xtick.labelsize"] = FONT_SIZE_SMALL
    mpl.rcParams["ytick.labelsize"] = FONT_SIZE_SMALL
    mpl.rcParams["legend.fontsize"] = FONT_SIZE_SMALL
    mpl.rcParams["lines.linewidth"] = 1.0
    mpl.rcParams["axes.linewidth"] = 0.8
    mpl.rcParams["grid.linewidth"] = 0.5
    mpl.rcParams["grid.alpha"] = 0.5
    mpl.rcParams["xtick.major.width"] = 0.8
    mpl.rcParams["ytick.major.width"] = 0.8

    # --- Legend Settings ---
    mpl.rcParams["legend.frameon"] = True
    mpl.rcParams["legend.framealpha"] = 0.75
    mpl.rcParams["legend.edgecolor"] = "black"  # Or 'inherit'
    mpl.rcParams["legend.fancybox"] = False

    # --- Axis Ticks ---
    mpl.rcParams["xtick.direction"] = "out"
    mpl.rcParams["ytick.direction"] = "out"
    mpl.rcParams["xtick.top"] = False
    mpl.rcParams["ytick.right"] = False
