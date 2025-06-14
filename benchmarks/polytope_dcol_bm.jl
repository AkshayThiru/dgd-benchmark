using Pkg
const ENV_PATH = joinpath(@__DIR__, "../extern/julia_env")
println(ENV_PATH)
Pkg.activate(ENV_PATH)

include("../benchmarks/helpers/random_utils.jl")
include("../benchmarks/helpers/benchmark_result.jl")
include("../benchmarks/helpers/mesh_loader.jl")

using Statistics
using BenchmarkTools
using ProgressBars


nstep = 10
npair = 10
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

rng = Random.MersenneTwister(5489)
nrows = [round(Int, 2^p) for p in range(3, 8.5, length=nstep)]
position_range = 5.0

nruns = nstep * npair * npose
res_arr = BenchmarkResultArray(nruns)
pbar = ProgressBar(total = nruns)
for nrow in nrows
  for i in 1:npair
    A₁, b₁ = generate_random_polyhedron(rng, nrow)
    A₂, b₂ = generate_random_polyhedron(rng, nrow)
    set1 = dc.Polytope(A₁, b₁)
    set2 = dc.Polytope(A₂, b₂)
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
      add_result!(res_arr, solve_time * 1e6, prim_dual_gap, prim_feas_err, nrow)
      ProgressBars.update(pbar)
    end
  end
end

save_to_feather_file(res_arr, joinpath(log_path, "polytope_bm__cold_dcol.feather"), true)
println("Avg. solve time (us): ", mean(res_arr.solve_times))
println("Max. prim dual gap  : ", maximum(res_arr.prim_dual_gaps))
println("Max. prim feas err  : ", maximum(res_arr.prim_feas_errs))
