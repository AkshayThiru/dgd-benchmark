import argparse
import os

import matplotlib.pyplot as plt

import plotting_tools as pt
from plotting_tools.common_plotting import (generate_combined_box_plots,
                                            generate_filled_line_plot,
                                            generate_row_box_plots, save_plot)
from plotting_tools.data_reader import (calculate_masked_statistics,
                                        read_feather_files_from_directory)
from plotting_tools.ieee_config import IEEE_SINGLE_COLUMN_WIDTH


def generate_latex_table_rows(
    data: list[list[dict]],
    exps: list[str],
    methods: list[str],
    fields: list[str],
    float_precision: int = 4,
) -> str:
    latex_rows = []
    latex_rows.append(" & ".join(["", ""] + methods) + r" \\")
    for e, row in enumerate(data):
        for field in fields:
            formatted_cells = []
            formatted_cells.append(exps[e].replace("\n", " ") + " & " + field)
            for cell in row:
                if not cell:
                    formatted_cells.append("-")
                else:
                    formatted_cells.append(
                        f"${cell['solve_time'][field]:.{float_precision}f}$"
                    )
            latex_rows.append(" & ".join(formatted_cells) + r" \\")
    return "\n".join(latex_rows)


def main():
    parser = argparse.ArgumentParser(description="Generate benchmark box plots.")
    parser.add_argument(
        "log_directory",
        type=str,
        help="Path to the directory containing Feather log files.",
    )
    args = parser.parse_args()
    log_dir = args.log_directory
    output_dir = "plots/output"
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    df_dicts = []
    exps = []
    methods = ["dgd", "dcf", "ie", "inc", "dcol"]
    method_labels = [r"Our", r"DCF", r"IE", r"Inc", r"DCol"]
    # Get data for combined box plot
    #   Primitive convex sets, cold-start
    df_dicts.append(
        read_feather_files_from_directory(log_dir, "primitive_bm__cold_", methods)
    )
    exps.append(r"Primitives," "\n" r"cold start")
    #   Primitive convex sets, warm-start
    df_dicts.append(
        read_feather_files_from_directory(log_dir, "primitive_bm__warm_", methods)
    )
    exps.append(r"Primitives," "\n" r"warm start")
    #   Mesh sets, cold-start
    df_dicts.append(
        read_feather_files_from_directory(log_dir, "mesh_bm__cold_", methods)
    )
    exps.append(r"Meshes," "\n" r"cold start")
    #   Mesh sets, warm-start
    df_dicts.append(
        read_feather_files_from_directory(log_dir, "mesh_bm__warm_", methods)
    )
    exps.append(r"Meshes," "\n" r"warm start")

    # Plots
    pt.set_ieee_single_column_plot_style()

    #   Solution time and iterations
    fig, axs = generate_combined_box_plots(
        # Data
        df_dicts,
        exps,
        methods,
        method_labels,
        ["solve_time", "iter"],
        [r"Solution time ($\mu s$)", r"Iterations"],
        # Figure
        IEEE_SINGLE_COLUMN_WIDTH,
        3.5,
        # Box properties
        whis=[0.01, 99.99],
    )
    for i in range(axs.shape[1]):
        axs[0, i].set_yscale("log")
        axs[1, i].set_yscale("log")
        axs[1, i].set_ylim([0.8, 200])
    axs[0, 1].margins(x=0.1, y=0.1)
    axs[1, 1].margins(x=0.1, y=0.1)

    plt.show()
    save_plot(fig, os.path.join(output_dir, "solution_time_iters.png"), dpi=750)

    #   DSF: solutions and iterations
    #   DSF sets, cold-start
    df_dict = read_feather_files_from_directory(log_dir, "dsf_bm__cold_", methods)
    fig, axs = generate_row_box_plots(
        # Data
        df_dict,
        methods,
        method_labels,
        ["solve_time", "iter"],
        [r"Solution time ($\mu s$)", r"Iterations"],
        # Figure
        IEEE_SINGLE_COLUMN_WIDTH,
        1.2,
        # Box properties
        box_width=1.0,
        gap=0.25,
        whis=[0.01, 99.99],
    )
    axs[0, 0].set_yscale("log")
    axs[0, 1].set_yscale("log")
    axs[0, 1].set_ylim([6, 110])

    plt.show()
    save_plot(fig, os.path.join(output_dir, "solution_time_iters_dsf.png"), dpi=750)

    df_dicts.append(df_dict)
    exps.append(r"DSF sets, cold start")

    #   Solution time vs polytope size
    df_dict = read_feather_files_from_directory(log_dir, "polytope_bm__cold_", ["dcol"])
    fig, ax = generate_filled_line_plot(
        # Data
        df_dict["dcol"],
        "polytope_size",
        "solve_time",
        r"Number of polytope hyperplanes",
        r"Solution time ($\mu s$)",
        # Figure
        IEEE_SINGLE_COLUMN_WIDTH,
        1.3,
        # Box properties
        percentiles=[0.1, 0.9],
    )
    ax.set_xscale("log")
    ax.set_yscale("log")

    plt.show()
    save_plot(fig, os.path.join(output_dir, "solution_time_polytope_size.png"), dpi=750)

    # Table of solution times
    # Create input for function
    data = []
    for e in range(len(exps)):
        row = []
        for method in methods:
            if method in df_dicts[e]:
                if method == "dcol":
                    df = calculate_masked_statistics(
                        df_dicts[e][method], ["solve_time"]
                    )
                else:
                    df = calculate_masked_statistics(
                        df_dicts[e][method], ["solve_time"], "optimal"
                    )
                row.append(df)
            else:
                row.append(dict())
        data.append(row)
    fields = ["p50", "p99.99"]
    print("\nLatex table:")
    print(generate_latex_table_rows(data, exps, methods, fields, float_precision=2))


if __name__ == "__main__":
    main()
