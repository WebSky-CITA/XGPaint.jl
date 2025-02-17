"""Thermal Sunyaev-Zel'dovich profiles, based on Battaglia et al. 2016. We also include a 
break model variant."""


struct Battaglia16ThermalSZProfile{T,C} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
end

struct BreakModel{T,C} <: AbstractGNFW{T}
    f_b::T
    cosmo::C
    alpha_break::T
    M_break::T
end

function Battaglia16ThermalSZProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    return Battaglia16ThermalSZProfile(f_b, cosmo)
end

function BreakModel(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, alpha_break::T=1.5, M_break::T=2.0*10^14) where {T <: Real}
    #alpha_break = 1.486 from Shivam P paper by Nate's sleuthing
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    return BreakModel(f_b, cosmo, alpha_break, M_break)
end

const ρ_crit_factor = uconvert(u"kg/m^3", 3u"km^2*Mpc^-2*s^-2" / (8π * constants.G))


function generalized_nfw(x, xc, α, β, γ)
    x̄ = x / xc
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end

function _generalized_scaled_nfw(x̄, α, β, γ)
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end


function get_params(::AbstractGNFW{T}, M_200, z) where T
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

function _nfw_profile_los_quadrature(x, xc, α, β, γ; zmax=1e5, rtol=eps(), order=9)
    x² = x^2
    scale = 1e9
    integral, err = quadgk(y -> scale * generalized_nfw(√(y^2 + x²), xc, α, β, γ),
                      0.0, zmax, rtol=rtol, order=order)
    return 2integral / scale
end

function dimensionless_P_profile_los(model::Battaglia16ThermalSZProfile{T}, r, M_200, z) where T
    par = get_params(model, M_200, z)
    R_200 = R_Δ(model, M_200, z, 200)
    x = r / angular_size(model, R_200, z)
    return par.P₀ * _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
end

function dimensionless_P_profile_los(model::BreakModel{T}, r, M_200, z) where T
    par = get_params(model, M_200, z)
    R_200 = R_Δ(model, M_200, z, 200)
    x = r / angular_size(model, R_200, z)
    if M_200 < model.M_break * M_sun
        return (par.P₀ * (M_200/(model.M_break*M_sun))^model.alpha_break * 
            _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ))
    else
        return par.P₀ * _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
    end
end

"""Line-of-sight integrated electron pressure"""
P_e_los(model, r, M_200c, z) = 0.5176 * P_th_los(model, r, M_200c, z)

"""Line-of-sight integrated thermal pressure"""
P_th_los(model, r, M_200c, z) = constants.G * M_200c * 200 * ρ_crit(model, z) * 
    model.f_b / 2 * dimensionless_P_profile_los(model, r, M_200c, z)

function compton_y(model, r, M_200c, z)
    return P_e_los(model, r, M_200c, z) * P_e_factor
end

(model::Battaglia16ThermalSZProfile)(r, M_200c, z) = compton_y(model, r, M_200c, z)
