# dgd-benchmark

Benchmarks for the [differentiable-growth-distance](https://github.com/HybridRobotics/differentiable-growth-distance) (DGD) library.


## Requirements

- C++17 compiler, CMake >= 3.15
- [Apache Arrow](https://arrow.apache.org/) C++ (for I/O)
- Python 3.10+ with:
  - pyarrow, matplotlib, seaborn, numpy, pandas
- Julia (for the .jl benchmarks)


## Setup

Install Python dependencies (for example, using a conda environment):
```sh
conda activate dgd-b
conda install -c conda-forge pyarrow matplotlib seaborn numpy pandas
# add project to PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```

Install Apache Arrow C++ (see link above). Ensure that Arrow's C++ headers & libs are discoverable by CMake.

Set up Julia environment (example):
```sh
# from repo root
julia scripts/setup_julia_env.jl
```

**Note:** For Mesh benchmarks, .obj files should be added to the `assets/` folder. For a mesh dataset, see the [YCB object dataset](https://www.ycbbenchmarks.com/).


## Build and run C++ benchmarks

From repository root, run:
```sh
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --parallel $(nproc)
```

Built benchmark binaries are located in `build/benchmarks/` (or the configured CMake output directory).

- Differentiable Support Function (DSF) set benchmark:
```sh
# usage: ./benchmarks/dsf_bm <assets_dir> <log_dir>
./benchmarks/dsf_bm /path/to/assets /path/to/logs
```

- Mesh benchmark:
```sh
# usage: ./benchmarks/mesh_bm <assets_dir> <log_dir>
./benchmarks/mesh_bm /path/to/assets /path/to/logs
```

- Primitive benchmark:
```sh
# usage: ./benchmarks/primitive_bm /path/to/logs
./benchmarks/primitive_bm /path/to/logs
```


## Run Julia benchmarks

From repository root, run:
```sh
# polytope (DCol) benchmark
julia benchmarks/polytope_dcol_bm.jl /path/to/logs

# primitive set benchmark
julia benchmarks/primitive_dcol_bm.jl /path/to/logs
```


## Plotting results

From the repository root, run:
```sh
# activate the conda environment
conda activate dgd-b
# ensure PYTHONPATH points to repository root so plotting utilities can import modules
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```

Use the Python plotting utilities in `plots/` to plot graphs:
```sh
python3 plots/bm_plot.py /path/to/logs
# output images are written to plots/output/
```

Additional flags for plotting:
- `--cp-lu`: Use LU decomposition results instead of Cramer's rule results for the cutting plane solver.
- `--warm-dual`: Use dual warm start results instead of primal warm start results.
- `--omit-trn`: Omit trust region Newton solver results.
