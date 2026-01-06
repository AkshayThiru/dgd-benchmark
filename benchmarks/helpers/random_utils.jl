using Random
using StaticArrays
using LinearAlgebra

import DifferentiableCollisions as dc


# --- Random Number Generators ---
@inline function random_range(rng::Random.AbstractRNG, min_val::Real, max_val::Real)
  return min_val + (max_val - min_val) * rand(rng)
end

@inline function random_position(rng::Random.AbstractRNG, abs_val::Real)
  return @SVector [random_range(rng, -abs_val, abs_val) for _ in 1:3]
end

@inline function random_unit_vector(rng::Random.AbstractRNG)
  n_raw = @SVector randn(rng, 3)
  return normalize(n_raw)
end

@inline function random_quaternion(rng::Random.AbstractRNG)
  q_raw = @SVector randn(rng, 4)
  return normalize(q_raw)
end


# --- Convex Primitive Parameter Ranges ---
struct CapsuleRanges{T}
  R::Tuple{T,T}
  L::Tuple{T,T}
end

struct CylinderRanges{T}
  R::Tuple{T,T}
  L::Tuple{T,T}
end

struct ConeRanges{T}
  H::Tuple{T,T}
  β::Tuple{T,T}
end

struct SphereRanges{T}
  R::Tuple{T,T}
end

struct EllipsoidRanges{T}
  L₁::Tuple{T,T}
  L₂::Tuple{T,T}
  L₃::Tuple{T,T}
end

struct PrimitiveParameterRanges{T}
  capsule::CapsuleRanges{T}
  cylinder::CylinderRanges{T}
  cone::ConeRanges{T}
  sphere::SphereRanges{T}
  ellipsoid::EllipsoidRanges{T}
end

function default_primitive_ranges(::Type{T}=Float64) where T
  PrimitiveParameterRanges{T}(
    CapsuleRanges{T}(
      (T(0.25 * 1e-2), T(0.25)), # R
      (T(0.5 * 1e-2), T(0.5)), # L
    ),
    CylinderRanges{T}(
      (T(0.25 * 1e-2), T(0.25)), # R
      (T(0.8 * 1e-2), T(0.8)), # L
    ),
    ConeRanges{T}(
      (T(0.5 * 1e-2), T(0.5)), # H
      (T(atan(0.5 * 1e-2)), T(atan(0.5 * 1e2))), # β
    ),
    SphereRanges{T}(
      (T(0.25 * 1e-2), T(0.25)) # R
    ),
    EllipsoidRanges{T}(
      (T(0.25 * 1e-2), T(0.25)), # L₁
      (T(0.25 * 1e-2), T(0.25)), # L₂
      (T(0.25 * 1e-2), T(0.25)), # L₃
    )
  )
end


# --- Convex Primitve Generator --
function randomize_parameters!(rng::Random.AbstractRNG, set::P, pr::PrimitiveParameterRanges{T}) where {P<:dc.AbstractPrimitive,T}
  error("No method implemented for randomizing $(typeof(set)) parameters.")
end

function randomize_parameters!(rng::Random.AbstractRNG, capsule::dc.Capsule{T}, pr::PrimitiveParameterRanges{T}) where T
  ranges = pr.capsule
  capsule.R = random_range(rng, ranges.R[1], ranges.R[2])
  capsule.L = random_range(rng, ranges.L[1], ranges.L[2])
  return capsule
end

function randomize_parameters!(rng::Random.AbstractRNG, cylinder::dc.Cylinder{T}, pr::PrimitiveParameterRanges{T}) where T
  ranges = pr.cylinder
  cylinder.R = random_range(rng, ranges.R[1], ranges.R[2])
  cylinder.L = random_range(rng, ranges.L[1], ranges.L[2])
  return cylinder
end

function randomize_parameters!(rng::Random.AbstractRNG, cone::dc.Cone{T}, pr::PrimitiveParameterRanges{T}) where T
  ranges = pr.cone
  cone.H = random_range(rng, ranges.H[1], ranges.H[2])
  cone.β = random_range(rng, ranges.β[1], ranges.β[2])
  return cone
end

function randomize_parameters!(rng::Random.AbstractRNG, sphere::dc.Sphere{T}, pr::PrimitiveParameterRanges{T}) where T
  ranges = pr.sphere
  sphere.R = random_range(rng, ranges.R[1], ranges.R[2])
  return sphere
end

function randomize_parameters!(rng::Random.AbstractRNG, ellipsoid::dc.Ellipsoid{T}, pr::PrimitiveParameterRanges{T}) where T
  ranges = pr.ellipsoid
  L₁ = random_range(rng, ranges.L₁[1], ranges.L₁[2])
  L₂ = random_range(rng, ranges.L₂[1], ranges.L₂[2])
  L₃ = random_range(rng, ranges.L₃[1], ranges.L₃[2])
  P = SMatrix{3,3}(Diagonal(SA[L₁, L₂, L₃]))
  ellipsoid.P = P
  ellipsoid.U = SMatrix{3,3}(cholesky(P).U)
  ellipsoid.F = eigen(P)
  return ellipsoid
end

function randomize_pose!(rng::Random.AbstractRNG, set::P, position_range::Real) where P<:dc.AbstractPrimitive
  set.r = random_position(rng, position_range)
  set.q = random_quaternion(rng)
end

function random_primitive_pair!(rng::Random.AbstractRNG, sets::Vector{<:dc.AbstractPrimitive}, pr::PrimitiveParameterRanges{T}) where T
  set1 = rand(rng, sets)
  set2 = rand(rng, sets)

  randomize_parameters!(rng, set1, pr)
  randomize_parameters!(rng, set2, pr)

  return (set1, set2)
end
