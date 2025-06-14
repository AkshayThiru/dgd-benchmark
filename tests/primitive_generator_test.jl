using Pkg
const ENV_PATH = joinpath(@__DIR__, "../extern/julia_env")
println(ENV_PATH)
Pkg.activate(ENV_PATH)

include("../benchmarks/helpers/random_utils.jl")


sets = dc.AbstractPrimitive[
  dc.Capsule(0.0, 0.0),
  dc.Cylinder(0.0, 0.0),
  dc.Cone(0.0, deg2rad(0.0)),
  dc.Sphere(0.0),
  dc.Ellipsoid(one(SMatrix{3,3,Float64,9})),
]

println("--- Initial state of sets ---")
for (i, set) in enumerate(sets)
  println("Set $i Type $(typeof(set))")
end

println("\n--- Randomizing pairs of existing primitives ---")
rng = Random.MersenneTwister(5489)
pr = default_primitive_ranges()
position_range = 5.0

for i in 1:10
  set1, set2 = random_primitive_pair!(rng, sets, pr)
  println("Pair $i: Types ($(typeof(set1)), $(typeof(set2)))")
  randomize_pose!(rng, set1, position_range)
end
