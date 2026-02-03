import argparse
import os

import plotting_tools as pt
from plotting_tools.common_plotting import (
    generate_combined_box_plots,
    generate_filled_line_plot,
    generate_line_plot_quantile,
    generate_row_box_plots,
    save_plot,
)
from plotting_tools.data_reader import (
    calculate_masked_statistics,
    read_feather_files_from_directory,
)
from plotting_tools.ieee_config import IEEE_SINGLE_COLUMN_WIDTH


def generate_latex_table_rows(
    data: list[list[dict]],
    exps: list[str],
    methods: list[str],
    fields: list[str],
    float_precision: int = 4,
) -> str:
    table_rows: list[list[str]] = []
    header: list[str] = ["", ""] + methods
    table_rows.append(header)

    for e, row_data in enumerate(data):
        exp_label = exps[e].replace("\n", " ")
        for field in fields:
            cells: list[str] = []
            cells.append(exp_label)
            cells.append(field)
            for cell in row_data:
                if not cell:
                    cells.append("-")
                else:
                    cells.append(f"${cell['solve_time'][field]:.{float_precision}f}$")
            table_rows.append(cells)

    # Compute max width per column
    ncols = 2 + len(methods)
    col_widths = [0] * ncols
    for r in table_rows:
        for c, val in enumerate(r):
            col_widths[c] = max(col_widths[c], len(val))

    # Format rows with padding: left-align first two columns, right-align numeric
    # columns
    lines: list[str] = []
    for row in table_rows:
        formatted_cells: list[str] = []
        for c, val in enumerate(row):
            if c < 2:
                formatted_cells.append(val.ljust(col_widths[c]))
            else:
                formatted_cells.append(val.rjust(col_widths[c]))
        lines.append(" & ".join(formatted_cells) + r" \\")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Generate benchmark box plots.")
    parser.add_argument(
        "log_directory",
        type=str,
        help="Path to the directory containing Feather log files.",
    )
    parser.add_argument(
        "--cp-lu",
        action="store_true",
        help="Cp solver: plot LU results, instead of Cramer",
    )
    parser.add_argument(
        "--warm-dual",
        action="store_true",
        help="Cp solver: plot dual warm start results, instead of primal warm start",
    )
    parser.add_argument(
        "--omit-trn", action="store_true", help="Omit Trn solver results"
    )
    args = parser.parse_args()
    log_dir = args.log_directory
    plot_cp_lu = args.cp_lu
    plot_warm_dual = args.warm_dual
    omit_trn = args.omit_trn
    output_dir = "plots/output"
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    df_dicts = []
    exps = []
    # Methods
    methods = []
    if plot_cp_lu:
        methods += ["dgd_cp_lu"]
        methods += ["dgd_dual_cp_lu"] if plot_warm_dual else ["dgd_primal_cp_lu"]
    else:
        methods += ["dgd_cp_cramer"]
        methods += (
            ["dgd_dual_cp_cramer"] if plot_warm_dual else ["dgd_primal_cp_cramer"]
        )
    methods += [] if omit_trn else ["dgd_trn"]
    methods += ["dcf", "ie", "inc", "dcol"]
    # Method labels
    # dgd_name = r"DGD"  # NOTE
    dgd_name = r"Our"  # For double blind review
    method_labels = []
    if omit_trn:
        method_labels += [dgd_name, dgd_name]  # Same name for cold and warm start
    else:
        method_labels += [
            dgd_name + "\n" r"(Cp)",
            dgd_name + "\n" r" (Cp)",
            dgd_name + "\n" r" (Trn)",
        ]
    method_labels += [r"DCF", r"IE", r"Inc", r"DCol"]

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
        IEEE_SINGLE_COLUMN_WIDTH,  # 1.5 * IEEE_SINGLE_COLUMN_WIDTH
        2.5,  # 3.5
        # Box properties
        whis=[0.01, 99.99],
    )
    for i in range(axs.shape[1]):
        axs[0, i].set_yscale("log")
        axs[1, i].set_yscale("log")
        axs[1, i].set_ylim([0.8, 200])
    axs[0, 1].margins(x=0.1, y=0.1)
    axs[1, 1].margins(x=0.1, y=0.1)

    # plt.show()  # NOTE
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

    # plt.show()  # NOTE
    save_plot(fig, os.path.join(output_dir, "solution_time_iters_dsf.png"), dpi=750)

    df_dicts.append(df_dict)
    exps.append(r"DSF sets, cold start")

    #   Solution time vs polytope size
    method_keys = ["dcol", "dgd_polytope", "dgd_mesh"]
    method_labels = [r"DCol", dgd_name, dgd_name + r" (hill-climbing)"]
    df_dict = read_feather_files_from_directory(
        log_dir, "polytope_bm__cold_", method_keys
    )
    fig, ax = generate_filled_line_plot(
        # Data
        df_dict,
        method_keys,
        method_labels,
        "polytope_size",
        "solve_time",
        r"Number of polytope vertices",
        r"Solution time ($\mu s$)",
        # Figure
        IEEE_SINGLE_COLUMN_WIDTH,
        1.75,
        # Colors
        line_colors=["#1f77b4", "#ff7f0e", "#2ca02c"],
        fill_colors=["#1f77b4", "#ff7f0e", "#2ca02c"],
        alpha=0.2,
        # Markers
        markersize=0.5,
        # Box properties
        percentiles=[0.0, 1.0],
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim([3e-1, 1e4])
    if len(method_keys) > 1:
        ax.legend(loc="upper right")

    # plt.show()  # NOTE
    save_plot(fig, os.path.join(output_dir, "solution_time_polytope_size.png"), dpi=750)

    #   Convergence error vs iteration
    method_keys = ["primitive", "mesh_1e-2"]
    method_suffixes = ["primitive", "mesh"]
    df_dict = read_feather_files_from_directory(
        log_dir, "convergence_bm__", method_keys
    )
    for i in range(len(method_keys)):
        fig, ax = generate_line_plot_quantile(
            # Data
            df_dict[method_keys[i]],
            r"Iteration",
            r"Relative gap",
            # Figure
            IEEE_SINGLE_COLUMN_WIDTH,
            1.3,
        )
        ax.set_yscale("log")
        ax.set_ylim([1e-9, 1e3])

        # plt.show()  # NOTE
        save_plot(
            fig,
            os.path.join(output_dir, "convergence_rate_" + method_suffixes[i] + ".png"),
            dpi=750,
        )

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
