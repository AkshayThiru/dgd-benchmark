using Arrow


mutable struct BenchmarkResultArray
  solve_times::Vector{Float64}
  prim_dual_gaps::Vector{Float64}
  prim_feas_errs::Vector{Float64}
  iters::Vector{Int}
  polytope_sizes::Vector{Int}
  capacity::Int32
  current_size::Int32
end

function BenchmarkResultArray(capacity::Int)
  solve_times = Vector{Float64}(undef, capacity)
  prim_dual_gaps = Vector{Float64}(undef, capacity)
  prim_feas_errs = Vector{Float64}(undef, capacity)
  iters = Vector{Int}(undef, capacity)
  polytope_sizes = Vector{Int}(undef, capacity)

  return BenchmarkResultArray(solve_times, prim_dual_gaps, prim_feas_errs, polytope_sizes, iters, capacity, 0)
end

function add_result!(res_arr::BenchmarkResultArray, solve_time::Float64, prim_dual_gap::Float64, prim_feas_err::Float64, iter::Int; polytope_size::Int = -1, )
  if res_arr.current_size < res_arr.capacity
    idx = res_arr.current_size + 1
    res_arr.solve_times[idx] = solve_time
    res_arr.prim_dual_gaps[idx] = prim_dual_gap
    res_arr.prim_feas_errs[idx] = prim_feas_err
    res_arr.iters[idx] = iter
    res_arr.polytope_sizes[idx] = polytope_size
    res_arr.current_size += 1
  end
end

function save_to_feather_file(res_arr::BenchmarkResultArray, filename::String, is_polytope_bm::Bool)
  if is_polytope_bm
    arrow_table = (
      solve_time = res_arr.solve_times,
      prim_dual_gap = res_arr.prim_dual_gaps,
      prim_feas_err = res_arr.prim_feas_errs,
      iter = res_arr.iters,
      polytope_size = res_arr.polytope_sizes,
    )
  else
    arrow_table = (
      solve_time = res_arr.solve_times,
      prim_dual_gap = res_arr.prim_dual_gaps,
      prim_feas_err = res_arr.prim_feas_errs,
      iter = res_arr.iters,
    )
  end

  try
    Arrow.write(filename, arrow_table)
    println("Successfully wrote data to '$filename'")
  catch e
    @error "Error writing Arrow file: " exception=(e, catch_backtrace())
  end
end
