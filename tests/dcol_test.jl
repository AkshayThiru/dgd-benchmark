using Pkg
const ENV_PATH = joinpath(@__DIR__, "../extern/julia_env")
println(ENV_PATH)
Pkg.activate(ENV_PATH)

import DifferentiableCollisions as dc
using Test
# using BenchmarkTools
using LinearAlgebra
using StaticArrays

import FiniteDiff
import Random
Random.seed!(1234)


C₁ = dc.Cone(2.0, deg2rad(22))
C₁.r = 0.3 * (@SVector randn(3))
C₁.q = normalize((@SVector randn(4)))

C₂ = dc.Capsule(0.3, 1.2)
C₂.r = (@SVector randn(3))
C₂.q = normalize((@SVector randn(4)))


α₁, x₁ = dc.proximity(C₁, C₂; verbose=false, pdip_tol=1e-6)
α₂, x₂ = dc.proximity(C₂, C₁; verbose=false, pdip_tol=1e-6)
@test abs(α₁ - α₂) < 1e-8
@test norm(x₁ - x₂) < 1e-8
# @btime dc.proximity(C₁, C₂)
#info α₁


α₁, dα₁_dθ = dc.proximity_gradient(C₁, C₂; pdip_tol=1e-6)
α₂, dα₂_dθ = dc.proximity_gradient(C₂, C₁; pdip_tol=1e-6)

@test norm([dα₂_dθ[8:14]; dα₂_dθ[1:7]] - dα₁_dθ) < 1e-4
# @btime dc.proximity_gradient(C₁, C₂)


α₁, x₁, ∂z₁_∂θ = dc.proximity_jacobian(C₁, C₂; pdip_tol=1e-6)
α₂, x₂, ∂z₂_∂θ = dc.proximity_jacobian_slow(C₁, C₂; pdip_tol=1e-6)
α₃, x₃, ∂z₃_∂θ = dc.proximity_jacobian(C₂, C₁; pdip_tol=1e-6)

@test norm(∂z₁_∂θ - ∂z₂_∂θ) < 1e-5
@test norm([∂z₃_∂θ[:, 8:14] ∂z₃_∂θ[:, 1:7]] - ∂z₁_∂θ) < 1e-3
# @btime dc.proximity_jacobian(C₁, C₂)

# Check derivatives.
function fd_α(C₁, C₂, θ₁, θ₂)
  C₁.r = θ₁[1:3]
  C₁.q = θ₁[4:7]
  C₂.r = θ₂[1:3]
  C₂.q = θ₂[4:7]
  α, x = dc.proximity(C₁, C₂; pdip_tol=1e-6)
  [x; α]
end

J = FiniteDiff.finite_difference_jacobian(θ -> fd_α(C₁, C₂, θ[1:7], θ[8:14]), [C₁.r; C₁.q; C₂.r; C₂.q])

@test norm(J - ∂z₁_∂θ) < 1e-2
@test norm(J - ∂z₂_∂θ) < 1e-2
@test norm(J[4, :] - dα₁_dθ) < 1e-3
