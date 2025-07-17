import os
from typing import Optional, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd
import seaborn as sns


def generate_box_plot(
    # Data
    df: pd.DataFrame,  # DataFrame with columns: x_column_name, y_column_name
    x_column_name: str,
    y_column_name: str,
    x_label: str,
    y_label: str,
    # Figure
    fig_width: float,  # in inches
    fig_height: float,  # in inches
    # Colors
    fill: bool = False,
    # Box properties
    box_width: float = 1.0,
    gap: float = 0.1,
    whis: Union[
        float, list[float]
    ] = 0.15,  # whiskers at whis * IQR or at [percentile_low, percentile_high]
    fliersize: Optional[float] = 0.2,  # Size of dots used for outliers
    box_linewidth: float = 1.0,
) -> tuple[mpl.figure.Figure, mpl.axes.Axes]:
    fig, ax = plt.subplots(
        figsize=(fig_width, fig_height), dpi=400, layout="constrained"
    )
    sns.boxplot(
        # Data
        x=x_column_name,
        y=y_column_name,
        # hue=exp_column_name,
        data=df,
        ax=ax,
        # Colors
        # palette=palette,
        # saturation=saturation,
        fill=fill,
        # Box properties
        width=box_width,
        gap=gap,
        whis=whis,
        fliersize=fliersize,
        linewidth=box_linewidth,
    )
    ax.set_xlabel(x_label)
    ax.tick_params(axis="x", labelbottom=True)
    ax.set_ylabel(y_label)
    ax.tick_params(axis="y", labelleft=True)
    # Other properties
    ax.margins(x=0.05, y=0.05)
    ax.grid(axis="y")
    return fig, ax


def generate_filled_line_plot(
    # Data
    df: pd.DataFrame,  # DataFrame with columns: x_column_name, y_column_name
    x_column_name: str,
    y_column_name: str,
    x_label: str,
    y_label: str,
    # Figure
    fig_width: float,  # in inches
    fig_height: float,  # in inches
    # Colors
    line_color: str = "blue",
    fill_color: str = "skyblue",
    # Line properties
    marker: str = "o",
    markersize: float = 1,
    # Fill properties
    percentiles: list[float] = [0.01, 0.99],
    alpha: float = 0.4,
) -> tuple[mpl.figure.Figure, mpl.axes.Axes]:
    summary_df = (
        df.groupby(x_column_name)[y_column_name]
        .agg(
            median=lambda x: x.quantile(0.5),
            plow=lambda x: x.quantile(percentiles[0]),
            phigh=lambda x: x.quantile(percentiles[1]),
        )
        .reset_index()
    )

    fig, ax = plt.subplots(
        figsize=(fig_width, fig_height), dpi=400, layout="constrained"
    )
    ax.plot(
        summary_df[x_column_name],
        summary_df["median"],
        marker=marker,
        markersize=markersize,
        ls="-",
        color=line_color,
    )
    ax.fill_between(
        summary_df[x_column_name],
        summary_df["plow"],
        summary_df["phigh"],
        color=fill_color,
        alpha=alpha,
    )

    ax.set_xlabel(x_label)
    ax.tick_params(axis="x", labelbottom=True)
    ax.set_ylabel(y_label)
    ax.tick_params(axis="y", labelleft=True)
    # Other properties
    ax.margins(x=0.05, y=0.05)
    ax.grid(axis="y")
    return fig, ax


def generate_combined_box_plots(
    # Data
    df_dicts: list[
        dict
    ],  # A list of dictionaries: [{'method_A': df_1, 'method_B': df_2}, {'method_A': df_3, 'method_B': df_4}]
    df_dict_labels: list[
        str
    ],  # A list of labels corresponding to each df_dict, e.g., ['Experiment X', 'Experiment Y']
    method_keys: list[str],  # e.g., ['method_A', 'method_B']
    method_labels: Optional[list[str]],  # Labels for the methods
    value_column_names: list[str],
    y_labels: list[str],
    # Figure
    fig_width: float,  # in inches
    fig_height: float,  # in inches
    # Color
    palette: str = "muted",
    saturation: float = 1.0,
    fill: bool = False,
    # Box properties
    box_width: float = 1.0,
    gap: float = 0.1,
    whis: Union[
        float, list[float]
    ] = 0.15,  # whiskers at whis * IQR or at [percentile_low, percentile_high]
    fliersize: Optional[float] = 0.2,  # Size of dots used for outliers
    box_linewidth: float = 0.75,
) -> tuple[mpl.figure.Figure, np.ndarray]:
    column_width_ratios = []
    for df_dict in df_dicts:
        column_width_ratios.append(len(df_dict))

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=400, layout="constrained")
    nrows = len(value_column_names)
    ncols = len(df_dicts)
    gs = fig.add_gridspec(
        nrows=nrows, ncols=ncols, width_ratios=column_width_ratios, wspace=0.1
    )
    axs = np.empty((nrows, ncols), dtype=object)

    for r_idx in range(nrows):
        value_column_name = value_column_names[r_idx]
        y_label = y_labels[r_idx]

        for c_idx, exp in enumerate(df_dict_labels):
            # Determine shared y-axis for the current row
            if c_idx == 0:
                axs[r_idx, c_idx] = fig.add_subplot(gs[r_idx, c_idx])
            else:
                axs[r_idx, c_idx] = fig.add_subplot(
                    gs[r_idx, c_idx], sharey=axs[r_idx, 0]
                )

            ax = axs[r_idx, c_idx]
            exp_method_data = []
            method_ordered = []
            for m, method in enumerate(method_keys):
                if method in df_dicts[c_idx]:
                    df = df_dicts[c_idx][method].copy()
                    df["__method__"] = method_labels[m]
                    df_subset = df[["__method__", value_column_name]]
                    exp_method_data.append(df_subset)
                    method_ordered.append(method_labels[m])
            exp_method_df = pd.concat(exp_method_data, ignore_index=True)
            # print(exp_method_df.info())

            palette_dict = dict(
                zip(method_labels, sns.color_palette(palette, len(method_labels)))
            )

            sns.boxplot(
                # Data
                x="__method__",
                y=value_column_name,
                # hue='__method__',
                data=exp_method_df,
                ax=ax,
                # Color
                # palette=palette_dict,
                # saturation=saturation,
                fill=fill,
                # Box properties
                dodge=True,
                order=method_ordered,
                width=box_width,
                gap=gap,
                whis=whis,
                fliersize=fliersize,
                linewidth=box_linewidth,
            )
            if r_idx == 0:
                ax.set_xlabel("")
                ax.tick_params(axis="x", labelbottom=False)
            else:
                ax.set_xlabel(exp)
                ax.tick_params(axis="x", labelbottom=True)
            if c_idx == 0:
                ax.set_ylabel(y_label)
            else:
                ax.set_ylabel("")
                ax.tick_params(axis="y", labelleft=False)
            # Other properties
            ax.margins(x=0.05, y=0.05)
            ax.grid(axis="y")

    return fig, axs


def generate_row_box_plots(
    # Data
    df_dict: dict,  # Dictionary: {'method_A': df_1, 'method_B': df_2}
    method_keys: list[str],  # e.g., ['method_A', 'method_B']
    method_labels: Optional[list[str]],  # Labels for the methods
    value_column_names: list[str],
    y_labels: list[str],
    # Figure
    fig_width: float,  # in inches
    fig_height: float,  # in inches
    # Color
    palette: str = "muted",
    saturation: float = 1.0,
    fill: bool = False,
    # Box properties
    box_width: float = 1.0,
    gap: float = 0.1,
    whis: Union[
        float, list[float]
    ] = 0.15,  # whiskers at whis * IQR or at [percentile_low, percentile_high]
    fliersize: Optional[float] = 0.2,  # Size of dots used for outliers
    box_linewidth: float = 0.75,
) -> tuple[mpl.figure.Figure, np.ndarray]:
    ncols = len(y_labels)
    column_width_ratios = [1.0 / ncols] * ncols

    fig = plt.figure(figsize=(fig_width, fig_height), dpi=400, layout="constrained")
    gs = fig.add_gridspec(
        nrows=1, ncols=ncols, width_ratios=column_width_ratios, wspace=0.1
    )
    axs = np.empty((1, ncols), dtype=object)

    for c_idx in range(ncols):
        value_column_name = value_column_names[c_idx]
        y_label = y_labels[c_idx]
        # Determine shared y-axis for the current row
        axs[0, c_idx] = fig.add_subplot(gs[0, c_idx])

        ax = axs[0, c_idx]
        exp_method_data = []
        method_ordered = []
        for m, method in enumerate(method_keys):
            if method in df_dict:
                df = df_dict[method].copy()
                df["__method__"] = method_labels[m]
                df_subset = df[["__method__", value_column_name]]
                exp_method_data.append(df_subset)
                method_ordered.append(method_labels[m])
        exp_method_df = pd.concat(exp_method_data, ignore_index=True)
        # print(exp_method_df.info())

        palette_dict = dict(
            zip(method_labels, sns.color_palette(palette, len(method_labels)))
        )

        sns.boxplot(
            # Data
            x="__method__",
            y=value_column_name,
            # hue='__method__',
            data=exp_method_df,
            ax=ax,
            # Color
            # palette=palette_dict,
            # saturation=saturation,
            fill=fill,
            # Box properties
            dodge=True,
            order=method_ordered,
            width=box_width,
            gap=gap,
            whis=whis,
            fliersize=fliersize,
            linewidth=box_linewidth,
        )
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelbottom=True)
        ax.set_ylabel(y_label)
        # Other properties
        ax.margins(x=0.05, y=0.05)
        ax.grid(axis="y")

    return fig, axs


def save_plot(fig: mpl.figure.Figure, filename: str, dpi: float = 750) -> None:
    if filename:
        dir = os.path.dirname(filename)
        if dir:
            if not os.path.exists(dir):
                os.makedirs(dir)
            fig.savefig(filename, dpi=dpi, bbox_inches="tight")
            print(f"Plot saved to: {filename}")
        else:
            print("Directory not found for the filename")
