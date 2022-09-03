
import PhysicalConstants.CODATA2018 as constants
using Unitful
const M_sun = 1.98847e30u"kg"
const P_e_factor = constants.σ_e / (constants.m_e * constants.c_0^2)
using Cosmology
using QuadGK
using XGPaint


##

abstract type AbstractProfile{T} end

struct BattagliaProfile{T,C} <: AbstractProfile{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
end

function BattagliaProfile(Omega_c::T=0.2589, Omega_b::T=0.0486) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(Float64, h=0.6774, OmegaM=OmegaM)
    return BattagliaProfile(f_b, cosmo)
end


function generalized_nfw(x, xc, α, β, γ)
    x̄ = x / xc
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end

"""Line-of-sight integrated electron pressure"""
P_e_los(𝕡, M_200, z, r) = 0.5176 * P_th_los(𝕡, M_200, z, r)

"""Line-of-sight integrated thermal pressure"""
P_th_los(𝕡, M_200, z, r) = constants.G * M_200 * 200 * ρ_crit(𝕡, z) * 
    𝕡.f_b / 2 * dimensionless_P_profile_los(𝕡, M_200, z, r)

function compton_y(𝕡, M_200, z, r)
    return P_e_los(𝕡, M_200, z, r) * P_e_factor
end

function get_params(::BattagliaProfile{T}, M_200, z) where T
	z₁ = z + 1
	m = M_200 / (1e14M_sun)
	P₀ = 18.1 * m^0.154 * z₁^-0.758
	xc = 0.497 * m^-0.00865 * z₁^0.731
	β = 4.35 * m^0.0393 * z₁^0.415
	α = 1
    γ = -0.3
    β = γ - α * β  # Sigurd's conversion from Battaglia to standard NFW
    return (xc=T(xc), α=T(α), β=T(β), γ=T(γ), P₀=T(P₀))
end

function dimensionless_P_profile_los(𝕡::BattagliaProfile{T}, M_200, z, r) where T
    par = get_params(𝕡, M_200, z)
    R_200 = R_Δ(𝕡, M_200, z, 200)
    x = r / angular_size(𝕡, R_200, z)
    return par.P₀ * _tsz_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
end

function ρ_crit(𝕡, z)
    H_z = H(𝕡.cosmo, z)
    return uconvert(u"kg/m^3", 3H_z^2 / (8π * constants.G))
end

function R_Δ(𝕡, M_Δ, z, Δ=200)
    return (M_Δ / (4π/3 * Δ * ρ_crit(𝕡, z)))^(1/3)
end

function angular_size(𝕡::AbstractProfile, physical_size, z)
    d_A = angular_diameter_dist(𝕡.cosmo, z)

    # convert both to the same units and strip units for atan
    phys_siz_unitless = ustrip(uconvert(unit(d_A), physical_size))
    d_A_unitless = ustrip(d_A)
    return atan(phys_siz_unitless, d_A_unitless)
end

_tsz_y₁(x, _a) = (x*(_a+1))^(1/(_a+1))
_tsz_x₁(y, _a) = y^(_a+1)/(_a+1)
function _tsz_profile_los_quadrature(x, xc, α, β, γ; zmax=1e5, _a=8, rtol=sqrt(eps()))
    x² = x^2
    integral, err = quadgk(y -> y^_a * generalized_nfw(√(_tsz_x₁(y,_a)^2 + x²), xc, α, β, γ),
                      0.0, _tsz_y₁(zmax,_a), rtol=rtol)
    return 2integral
end


##
using BenchmarkTools, Test
p = BattagliaProfile()
@test abs(1 - compton_y(p, 3e15M_sun, 1.0, 0.0) / 0.000791435242101093) < 1e-5
@test abs(1 - compton_y(p, 5e11M_sun, 0.5, 1e-5) / 2.038204776991399e-08) < 1e-5
@test abs(1 - compton_y(p, 5e11M_sun, 1.0, 2e-5) / 1.6761760790073166e-08) < 1e-5
@test abs(1 - compton_y(p, 1e12M_sun, 2.0, 3e-5) /  4.646234867480562e-08) < 1e-5

##

compton_y(p, 1e12M_sun, 2.0, 3e-5)

##
rr = LinRange(0., deg2rad(4/60), 100)
plot(rad2deg.(rr) * 60,
    [compton_y(p, 1e12M_sun, 2.0, r) for r in rr],
    ylim=(0., 1e-10))

##
using Pixell
rft = RadialFourierTransform(n=256, pad=128)
ells = rft.revl
lbeam = exp.(-ells .* (ells .+ 1) .* deg2rad(1/60)^2)
plot(ells, lbeam)

##
# rprofs = [compton_y(p, 1e12M_sun, 2.0, r) for r in rft.r]
rprofs = [compton_y(p, 1e12M_sun, 2.0, r) for r in rft.r]
lprofs = real2harm(rft, rprofs)
lprofs .*= lbeam
rprofs2 = harm2real(rft, reverse(lprofs))


##
plot(rft.r[128:(end-128)], rprofs[128:(end-128)])
plot!(rft.r[128:(end-128)], rprofs2[128:(end-128)], xlim=(0., 20e-5), ls=:dash)

# lprofs .*= exp.(-rft.l .* (rft.l .+ 1) .* deg2rad(1/60)^2)

##

##
using DataInterpolations
r, rprofs = Pixell.unpad(rft, rft.r, rprofs)

##
plot(r, rprofs, xscale=:log)

##
plot(rft.r[64:end], rprofs[64:end], xscale=:log10)

##


