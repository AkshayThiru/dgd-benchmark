# Build C++ benchmarks
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . --parallel $(nproc)
cd ..

# Run C++ benchmarks
cd build
./benchmarks/dsf_bm ../assets/ ../logs/
./benchmarks/mesh_bm ../assets/ ../logs/
./benchmarks/primitive_bm ../logs/
cd ..


# Setup Julia environment
julia scripts/setup_julia_env.jl

# Run Julia benchmarks
julia benchmarks/polytope_dcol_bm.jl logs
julia benchmarks/primitive_dcol_bm.jl logs  # NOTE


# Setup the Python environment
conda activate dgd
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Plot the results
python plots/bm_plot.py logs
