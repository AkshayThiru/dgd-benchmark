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
./benchmarks/polytope_dgd_bm ../logs/
cd ..


# Setup Julia environment
julia scripts/setup_julia_env.jl

# Run Julia benchmarks
julia benchmarks/polytope_dcol_bm.jl logs
julia benchmarks/primitive_dcol_bm.jl logs


# Build convergence benchmarks
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DDGD_EXTRACT_METRICS=ON ..
cmake --build . --parallel $(nproc)

# Run the convergence rate benchmarks
./benchmarks/convergence_bm ../logs/
cd ..


# Setup the Python environment
conda activate dgd
export PYTHONPATH="${PYTHONPATH}:$(pwd)"

# Plot the results
python plots/bm_plot.py logs
