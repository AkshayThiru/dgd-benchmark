using Pkg
const ENV_PATH = joinpath(@__DIR__, "../extern/julia_env")
println(ENV_PATH)
Pkg.activate(ENV_PATH)

include("../benchmarks/helpers/random_utils.jl")
include("../benchmarks/helpers/benchmark_result.jl")

using Statistics
using BenchmarkTools
using ProgressBars


npair = 100
npose = 100

pdip_tol = 1e-6
position_range = 5.0
bm_samples = 15
max_bm_time_per_run = 1 * 1e-3

if length(ARGS) < 1
  println("Usage: julia primitive_dcol_bm.jl <log_path>")
  exit(1)
end

log_path = ARGS[1]
if !isdir(log_path)
  println("Invalid log directory")
  exit(1)
end


sets = dc.AbstractPrimitive[
  dc.Capsule(0.0, 0.0),
  dc.Cylinder(0.0, 0.0),
  dc.Cone(0.0, deg2rad(0.0)),
  dc.Sphere(0.0),
  dc.Ellipsoid(one(SMatrix{3,3,Float64,9})),
]

rng = Random.MersenneTwister(5489)
pr = default_primitive_ranges()

nruns = npair * npose
res_arr = BenchmarkResultArray(nruns)
pbar = ProgressBar(total = nruns)
for i in 1:npair
  set1, set2 = random_primitive_pair!(rng, sets, pr)
  for j in 1:npose
    randomize_pose!(rng, set1, position_range)
    randomize_pose!(rng, set2, position_range)
    # Get problem matrices
    G_ort1, h_ort1, G_soc1, h_soc1 = dc.problem_matrices(set1, set1.r, set1.q)
    G_ort2, h_ort2, G_soc2, h_soc2 = dc.problem_matrices(set2, set2.r, set2.q)
    # Create and solve SOCP
    c, G, h, idx_ort, idx_soc1, idx_soc2 = dc.combine_problem_matrices(G_ort1, h_ort1, G_soc1, h_soc1,G_ort2, h_ort2, G_soc2, h_soc2)
    x, s, z = dc.solve_socp(c, G, h, idx_ort, idx_soc1, idx_soc2; verbose = false, pdip_tol = pdip_tol)
    # Compute solution errors
    prim_dual_gap = dot(s, z)
    prim_feas_err = norm(G*x + s - h)
    solve_time = @belapsed dc.proximity($set1, $set2; verbose = false, pdip_tol = pdip_tol) samples=bm_samples seconds=max_bm_time_per_run
    add_result!(res_arr, solve_time * 1e6, prim_dual_gap, prim_feas_err)
    update(pbar)
  end
end

save_to_feather_file(res_arr, joinpath(log_path, "primitive_bm__cold_dcol.feather"), false)
println("Avg. solve time (us): ", mean(res_arr.solve_times))
println("Max. prim dual gap  : ", maximum(res_arr.prim_dual_gaps))
println("Max. prim feas err  : ", maximum(res_arr.prim_feas_errs))
