import argparse

from plotting_tools.data_reader import (
    read_feather_files_from_directory,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Read Feather files.")
    parser.add_argument(
        "log_directory",
        type=str,
        help="Path to the directory containing Feather log files.",
    )
    parser.add_argument(
        "filenames",
        type=str,
        nargs="+",
        help="Path to the Feather file to read.",
    )
    args = parser.parse_args()
    log_dir = args.log_directory
    filenames = args.filenames

    df_dict = read_feather_files_from_directory(log_dir, "", filenames)
    for name, df in df_dict.items():
        print(f"Columns in DataFrame '{name}': {df.columns.tolist()}")


if __name__ == "__main__":
    main()
