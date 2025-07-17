import os
from typing import Optional

import pandas as pd


def read_feather_files_from_directory(
    dir_path: str, prefix: str, names: list[str]
) -> dict:
    if not os.path.isdir(dir_path):
        print(f"Error: Directory '{dir_path}' does not exist or is not a directory.")
        return {}

    df_dict = {}
    for name in names:
        filename = os.path.join(dir_path, prefix + name + ".feather")

        if not os.path.exists(filename):
            print(f"Warning: File '{filename}' not found. Skipping.")
            continue

        try:
            df = pd.read_feather(filename)
            df_dict[name] = df
            # print(f"Successfully read: '{filename}' (Shape: {df.shape})")
        except Exception as e:
            print(f"Error reading {filename}: {e}")
    return df_dict


def calculate_masked_statistics(
    df: pd.DataFrame, columns: list[str], mask_col: Optional[str] = None
) -> dict:
    df_processed = df.copy()
    if mask_col:
        if mask_col in df_processed.columns and pd.api.types.is_bool_dtype(
            df_processed[mask_col]
        ):
            df_processed = df_processed[df_processed[mask_col]]
        elif mask_col not in df_processed.columns:
            print(
                f"Warning: Masking column '{mask_col}' not found in DataFrame. Masking skipped."
            )
        else:
            print(
                f"Warning: Masking column '{mask_col}' is not of boolean type. Masking skipped."
            )

    # --- Calculate Numeric Statistics ---
    stats_data = {}
    for col in columns:
        if col in df_processed.columns and pd.api.types.is_numeric_dtype(
            df_processed[col]
        ):
            series = df_processed[col]
            stats_data[col] = {
                "mean": series.mean(),
                "std": series.std(),
                "p50": series.quantile(0.50),
                "p90": series.quantile(0.90),
                "p99.99": series.quantile(0.9999),
                "min": series.min(),
                "max": series.max(),
            }
        elif col not in df_processed.columns:
            print(f"Warning: Column '{col}' not found in DataFrame.")
        else:
            print(f"Warning: Column '{col}' is not numeric and will be skipped.")

    return stats_data


def generate_comparison_box_plot_df(
    df_dicts: list[
        dict
    ],  # A list of dictionaries: [{'method_A': df_1, 'method_B': df_2}, {'method_A': df_3, 'method_B': df_4}]
    df_dict_labels: list[
        str
    ],  # A list of labels corresponding to each df_dict, e.g., ['Experiment X', 'Experiment Y']
    method_keys: list[str],  # e.g., ['method_A', 'method_B']
    value_column_name: str,  # The name of the column within each DataFrame that holds the numerical values for plotting
    method_labels: Optional[list[str]] = None,  # Labels for the methods
) -> pd.DataFrame:
    if not df_dicts:
        print("Error: 'df_dicts' list is empty. No plots to generate.")
        return
    if not df_dict_labels or len(df_dicts) != len(df_dict_labels):
        print(
            "Error: 'df_dict_labels' must be provided and match the number of 'df_dicts'."
        )
        return
    if not method_keys or (
        method_labels is not None and len(method_keys) != len(method_labels)
    ):
        print(
            "Error: 'method_keys' list must be nonempty and must match the number of 'method_labels' (if provided)."
        )
        return

    # --- Prepare the data for plotting ---
    # We need to combine data from all relevant DataFrames into a single, long-form DataFrame
    # suitable for seaborn.

    combined_data = []
    for i, df_dict in enumerate(df_dicts):
        # For each experiment...
        exp = df_dict_labels[i]
        for m, method in enumerate(method_keys):
            # For each method...
            if method in df_dict and isinstance(df_dict[method], pd.DataFrame):
                df = df_dict[method].copy()
                if value_column_name not in df.columns:
                    print(
                        f"Warning: Column '{value_column_name}' not found in DataFrame for method '{method}' from experiment '{exp}'. Skipping this data slice."
                    )
                    continue

                # Add columns to identify the experiment and the method
                df["__exp__"] = exp
                df["__method__"] = (
                    method_labels[m] if method_labels is not None else method
                )
                df_subset = df[["__exp__", "__method__", value_column_name]]
                combined_data.append(df_subset)
            else:
                print(
                    f"Warning: Key '{method}' not found or is not a DataFrame in the data dictionary for label '{exp}'. Skipping."
                )

    if not combined_data:
        print(
            "No valid data found for plotting after combining. Check keys, column names, and dataframes."
        )
        return

    # We will use '__method__' on the x-axis, 'value_column_name' on the y-axis,
    # and '__exp__' for coloring/grouping the box plots.
    df_for_plotting = pd.concat(combined_data, ignore_index=True)

    # print(f"Data prepared for plotting. Combined DataFrame shape: {df_for_plotting.shape}")
    # print("First 5 rows of combined data:")
    # print(df_for_plotting.head())

    return df_for_plotting
