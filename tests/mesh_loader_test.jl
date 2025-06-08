using Pkg
const ENV_PATH = joinpath(@__DIR__, "../extern/julia_env")
println(ENV_PATH)
Pkg.activate(ENV_PATH)

include("../benchmarks/helpers/mesh_loader.jl")


const ASSETS_PATH = joinpath(@__DIR__, "../assets")

points = read_obj_vertices(joinpath(ASSETS_PATH, "rock_lowpoly.obj"))
A, b = compute_convex_hull_hrep(points)
println("Size of A (mesh): ", size(A))
println("Size of b (mesh): ", size(b))

A, b = generate_random_polyhedron(10)
println("Size of A (gen): ", size(A))
println("Size of b (gen): ", size(b))
