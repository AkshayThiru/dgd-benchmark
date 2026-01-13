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
  points = [SVector{3,Float64}(v[1], v[2], v[3]) for v in vertices]
  return points
end

function compute_convex_hull_hrep(points::Vector{SVector{3,Float64}})
  # Create V-representation polyhedron.
  V_rep = vrep(points)
  P = polyhedron(V_rep, CDDLib.Library())
  removevredundancy!(P)

  # Extract A and b where {x: a_i' * x <= b_i} is the polyhedron.
  A_rows = []
  b_vals = []
  for hspace in halfspaces(P)
    push!(A_rows, SVector{3,Float64}(hspace.a))
    push!(b_vals, hspace.β)
  end
  A = reduce(vcat, transpose.(A_rows))
  b = SVector{length(b_vals),Float64}(b_vals...)

  return A, b
end

function get_half_lengths(rng::Random.AbstractRNG, skew::Real)
  if !(1e-3 <= skew <= 1.0)
    error("skew must be in the range [1e-3, 1.0]")
  end

  log_skew = log(skew)
  rand_skew = exp(log_skew + rand(rng) * (0.0 - log_skew))
  half_lengths = SVector{3,Float64}(skew, rand_skew, 1.0)
  return half_lengths
end

function generate_random_polyhedron(rng::Random.AbstractRNG, nvert::Int, radius::Real=0.25)
  # Generate nvert random points on the ellipsoid surface.
  points = Vector{SVector{3,Float64}}(undef, nvert)
  for i in 1:nvert
    u = randn(rng, 3)
    u_normalized = radius * normalize(u)
    points[i] = SVector{3,Float64}(u_normalized)
  end

  # Set the center point as the origin.
  center = sum(points) / nvert
  points = [p - center for p in points]

  A, b = compute_convex_hull_hrep(points)
  nrows = size(A, 1)
  return SMatrix{nrows,3,Float64}(A), b
end
