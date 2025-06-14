# dgd-benchmark/scripts/setup_julia_env.jl

using Pkg

# This assumes the script is run from dgd-benchmark/scripts.
const EXTERN_DIR = joinpath(@__DIR__, "../extern")
const JULIA_ENV_DIR = joinpath(EXTERN_DIR, "julia_env")
const JULIA_DCOL_PATH = joinpath(EXTERN_DIR, "DifferentiableCollisions")

println("Setting up Julia environment in: $JULIA_ENV_DIR")
println("Using Julia repository from: $JULIA_DCOL_PATH")

# Activate the environment in julia_env.
Pkg.activate(JULIA_ENV_DIR)

# Add the external Julia repository as a development package.
Pkg.develop(path=JULIA_DCOL_PATH)

# Manually add core dependencies
Pkg.add("ForwardDiff")
Pkg.add("LinearAlgebra")
Pkg.add("MeshCat")
Pkg.add("Polyhedra")
Pkg.add("Printf")
Pkg.add("StaticArrays")
# Test dependencies
Pkg.add("BenchmarkTools")
Pkg.add("FiniteDiff")
# Mesh loader dependencies
Pkg.add("FileIO")
Pkg.add("MeshIO")
Pkg.add("GeometryBasics")
Pkg.add("Polyhedra")
Pkg.add("CDDLib")
# Logging dependencies
Pkg.add("Arrow")
# Others
Pkg.add("ProgressBars")

# Instantiate/resolve dependencies
Pkg.instantiate()

println("Julia environment setup complete. Manifest.toml generated.")
