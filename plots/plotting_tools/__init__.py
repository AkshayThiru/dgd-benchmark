from .ieee_config import (
    IEEE_DOUBLE_COLUMN_WIDTH,
    IEEE_SINGLE_COLUMN_WIDTH,
    apply_ieee_mpl_settings,
)


def set_ieee_single_column_plot_style():
    apply_ieee_mpl_settings(column_width="single")


def set_ieee_double_column_plot_style():
    apply_ieee_mpl_settings(column_width="double")


__all__ = [
    "set_ieee_single_column_plot_style",
    "set_ieee_double_column_plot_style",
    "IEEE_SINGLE_COLUMN_WIDTH",
    "IEEE_DOUBLE_COLUMN_WIDTH",
]
