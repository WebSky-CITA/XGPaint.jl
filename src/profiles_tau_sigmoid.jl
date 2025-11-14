"""Sigmoid-suppressed Battaglia tau (optical depth) profiles.
These profiles include a sigmoid suppression factor that smoothly transitions
from 0.99 at 4 R_200 to 0.01 at 6 R_200.
"""

# can be used with either angular (default) or physical units
#abstract type AbstractBattagliaTauProfile{T} <: AbstractGNFW{T} end

# -------------------------------------------------------------------------
# Sigmoid-suppressed tau models
# -------------------------------------------------------------------------

struct SigmoidBattagliaTauProfile{T,C} <: AbstractBattagliaTauProfile{T}
    f_b::T
    cosmo::C
end

struct SigmoidBattagliaTauProfilePhysical{T,C} <: AbstractBattagliaTauProfile{T}
    f_b::T
    cosmo::C
end

function SigmoidBattagliaTauProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774) where {T <: Real}
    OmegaM = Omega_b + Omega_c
    f_b    = Omega_b / OmegaM
    cosmo  = get_cosmology(T, h=h, Neff=3.046, OmegaM=OmegaM)
    return SigmoidBattagliaTauProfile{T, typeof(cosmo)}(f_b, cosmo)
end

function SigmoidBattagliaTauProfilePhysical(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774) where {T <: Real}
    OmegaM = Omega_b + Omega_c
    f_b    = Omega_b / OmegaM
    cosmo  = get_cosmology(T, h=h, Neff=3.046, OmegaM=OmegaM)
    return SigmoidBattagliaTauProfilePhysical{T, typeof(cosmo)}(f_b, cosmo)
end

# if angular, return the R200 size in radians
function object_size(model::SigmoidBattagliaTauProfile{T,C}, physical_size, z) where {T,C}
    d_A = angular_diameter_dist(model.cosmo, z)
    phys_siz_unitless = T(ustrip(uconvert(unit(d_A), physical_size)))
    d_A_unitless = T(ustrip(d_A))
    return atan(phys_siz_unitless, d_A_unitless)
end

# if physical, return the R200 size in Mpc
function object_size(::SigmoidBattagliaTauProfilePhysical{T,C}, physical_size, z) where {T,C}
    return physical_size
end



# -------------------------------------------------------------------------
# Sigmoid definition
# -------------------------------------------------------------------------
"""
    sigmoid_suppression(r)

Sigmoid suppression factor that equals 0.99 at r=4 and 0.01 at r=6.
Uses the form: f(r) = 1 / (1 + exp(k * (r - r0)))

Where k and r0 are solved from the boundary conditions:
- f(4) = 0.99 => 1/(1 + exp(k*(4-r0))) = 0.99
- f(6) = 0.01 => 1/(1 + exp(k*(6-r0))) = 0.01

This gives us:
- exp(k*(4-r0)) = 1/99 ≈ 0.0101
- exp(k*(6-r0)) = 99
- Taking the ratio: exp(2k) = 99^2 = 9801
- So k = ln(9801)/2 ≈ 4.599
- And r0 = 4 + ln(99)/k = 4 + ln(99)/4.599 ≈ 5
"""
function sigmoid_suppression(r::T) where {T <: Real}
    k  = T(2.3)
    r0 = T(5.0)
    return 1 / (1 + exp(k * (r - r0)))
end

# -------------------------------------------------------------------------
# Suppressed line-of-sight integration
# -------------------------------------------------------------------------
function _sigmoid_nfw_profile_los_quadrature(x, xc, α, β, γ; zmax=1e5, rtol=eps(), order=9)
    x² = x^2
    scale = 1e9
    integral, err = quadgk(
        y -> scale * sigmoid_suppression(√(y^2 + x²)) *
             XGPaint.generalized_nfw(√(y^2 + x²), xc, α, β, γ),
        0.0, zmax, rtol=rtol, order=order)
    return 2integral / scale
end

# -------------------------------------------------------------------------
# Density and number profiles with sigmoid suppression
# -------------------------------------------------------------------------
function rho_2d(model::SigmoidBattagliaTauProfile, r, m200c, z)
    par     = get_params(model, m200c, z)
    r200c   = R_Δ(model, m200c, z, 200)
    X       = r / object_size(model, r200c, z)
    rho_crit = ρ_crit_comoving_h⁻²(model, z)
    result  = par.P₀ * _sigmoid_nfw_profile_los_quadrature(X, par.xc, par.α, par.β, par.γ)
    return result * rho_crit * (r200c * (1 + z))
end

function rho_2d(model::SigmoidBattagliaTauProfilePhysical, r, m200c, z)
    par     = get_params(model, m200c, z)
    r200c   = R_Δ(model, m200c, z, 200)
    X       = r / object_size(model, r200c, z)
    rho_crit = ρ_crit_comoving_h⁻²(model, z)
    result  = par.P₀ * _sigmoid_nfw_profile_los_quadrature(X, par.xc, par.α, par.β, par.γ)
    return result * rho_crit * (r200c * (1 + z))
end
"""
function ne2d(model::SigmoidBattagliaTauProfile, r, m200c, z)
    me  = constants.ElectronMass
    mH  = constants.ProtonMass
    mHe = 4mH
    xH  = 0.76
    nH_ne  = 2xH / (xH + 1)
    nHe_ne = (1 - xH)/(2 * (1 + xH))
    factor = (me + nH_ne*mH + nHe_ne*mHe) / model.cosmo.h^2
    result = rho_2d(model, r, m200c, z)
    return result / factor
end

function ne2d(model::SigmoidBattagliaTauProfilePhysical, r, m200c, z)
    me  = constants.ElectronMass
    mH  = constants.ProtonMass
    mHe = 4mH
    xH  = 0.76
    nH_ne  = 2xH / (xH + 1)
    nHe_ne = (1 - xH)/(2 * (1 + xH))
    factor = (me + nH_ne*mH + nHe_ne*mHe) / model.cosmo.h^2
    result = rho_2d(model, r, m200c, z)
    return result / factor
end
"""
# -------------------------------------------------------------------------
# Compute τ (optical depth)
# -------------------------------------------------------------------------
"""
function compute_tau(model::SigmoidBattagliaTauProfile, r, m200c, z)
    return constants.ThomsonCrossSection * ne2d(model, r, m200c, z) + 0
end

function compute_tau(model::SigmoidBattagliaTauProfilePhysical, r, m200c, z)
    return constants.ThomsonCrossSection * ne2d(model, r, m200c, z) + 0
end
"""
# -------------------------------------------------------------------------
# Direct call overloads
# -------------------------------------------------------------------------
#(model::SigmoidBattagliaTauProfile)(r, m200c, z) =
#    compute_tau(model, r, m200c * M_sun, z)

#(model::SigmoidBattagliaTauProfilePhysical)(r, m200c, z) =
#    compute_tau(model, r, m200c * M_sun, z)
