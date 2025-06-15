using StaticArrays
using LinearAlgebra

import DifferentiableCollisions as dc


function custom_solve_socp(c::SVector{nx,T},
  G::SMatrix{ns,nx,T,nsnx},
  h::SVector{ns,T},
  idx_ort::SVector{n_ort,Ti},
  idx_soc1::SVector{n_soc1,Ti},
  idx_soc2::SVector{n_soc2,Ti};
  pdip_tol::T=1e-4
) where {nx,ns,nsnx,n_ort,n_soc1,n_soc2,T,Ti}
  x,s,z = dc.initialize(c,G,h,idx_ort,idx_soc1,idx_soc2)

  e = dc.gen_e(idx_ort, idx_soc1,idx_soc2)

  cone_degree = n_ort
  if n_soc1 > 0
    cone_degree += 1
  end
  if n_soc2 > 0
    cone_degree += 1
  end

  for main_iter = 1:50

    W = dc.calc_NT_scalings(s,z,idx_ort,idx_soc1,idx_soc2)

    # cache normalized variables
    λ = W*z
    λλ = dc.cone_product(λ,λ,idx_ort,idx_soc1,idx_soc2)

    # evaluate residuals
    rx = G'*z + c
    rz = s + G*x - h
    μ = dot(s,z)/(cone_degree)
    if μ < pdip_tol
      return x,s,z, main_iter
    end

    # affine step
    bx = -rx
    λ_ds = dc.inverse_cone_product(λ,-λλ,idx_ort,idx_soc1,idx_soc2)
    bz̃ = W\(-rz - W*(λ_ds))
    G̃ = W\G
    F = cholesky(Symmetric(G̃'*G̃))
    Δxa = F\(bx + G̃'*bz̃)
    Δza = W\(G̃*Δxa - bz̃)
    Δsa = W*(λ_ds - W*Δza)

    # linesearch on affine step
    αa = min(dc.linesearch(s,Δsa,idx_ort, idx_soc1,idx_soc2), dc.linesearch(z,Δza,idx_ort,idx_soc1,idx_soc2))
    ρ = dot(s + αa*Δsa, z + αa*Δza)/dot(s,z)
    σ = max(0, min(1,ρ))^3

    # centering and correcting step
    ds = -λλ - dc.cone_product(W\Δsa, W*Δza,idx_ort,idx_soc1,idx_soc2) + σ*μ*e
    λ_ds = dc.inverse_cone_product(λ,ds,idx_ort,idx_soc1,idx_soc2)
    bz̃ = W\(-rz - W*(λ_ds))
    Δx = F\(bx + G̃'*bz̃)
    Δz = W\(G̃*Δx - bz̃)
    Δs = W*(λ_ds - W*Δz)

    # final line search (.99 to avoid hitting edge of cone)
    α = min(1,0.99*min(dc.linesearch(s,Δs,idx_ort,idx_soc1,idx_soc2), dc.linesearch(z,Δz,idx_ort,idx_soc1,idx_soc2)))

    # take step
    x += α*Δx
    s += α*Δs
    z += α*Δz

  end
  error("pdip failed")

end
