using MeshIO
using FileIO
using Polyhedra
using StaticArrays
using LinearAlgebra
using GeometryBasics
using CDDLib
using Random

function read_obj_vertices(filepath::String)
  mesh = load(filepath)
  # Extract the coordinates of each point.
  vertices = GeometryBasics.coordinates(mesh)
  # Convert to a list of StaticArrays for convenience later.
  points = [SVector{3, Float64}(v[1], v[2], v[3]) for v in vertices]
  return points
end

function compute_convex_hull_hrep(points::Vector{SVector{3, Float64}})
  # Create V-representation polyhedron.
  V_rep = vrep(points)
  P = polyhedron(V_rep, CDDLib.Library())
  removevredundancy!(P)

  # Extract A and b where {x: a_i' * x <= b_i} is the polyhedron.
  A_rows = []
  b_vals = []
  for hspace in halfspaces(P)
    push!(A_rows, SVector{3, Float64}(hspace.a))
    push!(b_vals, hspace.β)
  end
  A = reduce(vcat, transpose.(A_rows))
  b = SVector{length(b_vals), Float64}(b_vals...)

  return A, b
end

function generate_random_polyhedron(nrows::Int, rng_seed::Int = 1234,
                                    max_norm_A::Real = 1.0,
                                    box_bound::Real = 5.0)
  Random.seed!(rng_seed)
  if nrows < 6
    error("To guarantee boundedness, nrows must be at least 6.")
  end

  A_rows = []
  b_vec = []

  for j in 1:3
    # Positive bound: x_j <= box_bound.
    ej_pos = zeros(3)
    ej_pos[j] = 1.0
    push!(A_rows, SVector{3}(ej_pos))
    push!(b_vec, box_bound)
    # Negative bound -x_j <= boux_bound.
    ej_neg = zeros(3)
    ej_neg[j] = -1.0
    push!(A_rows, SVector{3}(ej_neg))
    push!(b_vec, box_bound)
  end

  for i in 1:(nrows - 6)
    n = normalize(randn(3)) * rand() * max_norm_A
    b_i = rand() * (sqrt(3) * box_bound)

    push!(A_rows, SVector{3}(n))
    push!(b_vec, b_i)
  end

  A = SMatrix{nrows, 3, Float64}(reduce(vcat, transpose.(A_rows)))
  b = SVector{nrows, Float64}(b_vec)
  return A, b
end
