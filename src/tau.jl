
struct BattagliaTauProfile{T,C} <: AbstractProfile{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
end

function BattagliaTauProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, Neff=3.046, OmegaM=OmegaM)
    return BattagliaTauProfile(f_b, cosmo)
end

function get_params(::BattagliaTauProfile{T}, M_200, z) where T
	z₁ = z + 1
	m = M_200 / (1e14M_sun)
    P₀ = 4.e3 * m^0.29 * z₁^(-0.66)
	α =  0.88 * m^(-0.03) * z₁^0.19
	β = -3.83 * m^0.04 * z₁^(-0.025)
	xc = 0.5
    γ = -0.2
    return (xc=T(xc), α=T(α), β=T(β), γ=T(γ), P₀=T(P₀))
end



function ρ_crit_comoving_h⁻²(p, z)
    return  (ρ_crit(p, z) ) / (1+z)^3 / p.cosmo.h^2
end

function r200c_comoving(p, m200c, z)
    rho_crit = ρ_crit(p, z) / (1+z)^3 
    return cbrt(m200c/ (4π/3 * rho_crit * 200)) 
end


# returns a density, which we can check against Msun/Mpc² 
function rho_2d(p::BattagliaTauProfile, R_comoving_projected, m200c, z)
    par = get_params(p, m200c, z)
    r200c = r200c_comoving(p, m200c, z)
    X = R_comoving_projected / r200c
    rho_crit = ρ_crit_comoving_h⁻²(p, z)
<<<<<<< HEAD
    result = par.P₀ * XGPaint._nfw_profile_los_quadrature(X, par.xc, par.α, par.β, par.γ)
=======
    result = par.P₀ * XGPaint._tsz_profile_los_quadrature(X, par.xc, par.α, par.β, par.γ)
>>>>>>> main

    return result * rho_crit * r200c
end

function ne2d(p::BattagliaTauProfile, R_comoving_projected, m200c, z)
    me = constants.ElectronMass
    mH = constants.ProtonMass
    mHe = 4mH
    xH = 0.76
    nH_ne = 2xH / (xH + 1)
    nHe_ne = (1 - xH)/(2 * (1 + xH))
    factor = (me + nH_ne*mH + nHe_ne*mHe) / p.cosmo.h^2

    result = rho_2d(p, R_comoving_projected, m200c, z)  # (Msun/h) / (Mpc/h)^2
    return result / factor
end
