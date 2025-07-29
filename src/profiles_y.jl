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

struct Arnauld10ThermalSZProfile{T,C} <: AbstractGNFW{T}
    B::T
    cosmo::C
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

function Arnauld10ThermalSZProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, B::T=1.0) where {T <: Real}
    # Compute total matter density
    OmegaM = Omega_b + Omega_c
    # Obtain cosmology object which contains H(z) etc.
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    return Arnauld10ThermalSZProfile(B, cosmo)
end

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

# Overload get_params for the A10 model (adjust scaling as needed)
function get_params(::Arnauld10ThermalSZProfile{T}, M_200, z) where T
    P₀ = 8.130
    xc = 1 / 1.156
    α = 1.0620
    γ = -0.3292
    β = -5.4807
    return (xc=T(xc), α=T(α), β=T(β), γ=T(γ), P₀=T(P₀))
end


function _nfw_profile_los_quadrature(x, xc, α, β, γ; zmax=1e5, rtol=eps(), order=9)
    x² = x^2
    scale = 1e9
    integral, err = quadgk(y -> scale * generalized_nfw(√(y^2 + x²), xc, α, β, γ),
                      0.0, zmax, rtol=rtol, order=order)
    return 2integral / scale
end

#Check this normalisation
function A10_normalization(model::Arnauld10ThermalSZProfile{T}, m , z; B::T=one(T)) where {T <: Real}
    m_tilde = m / M_sun / B * model.cosmo.h # Convert mass to tilde-mass M_sun/h
    # println("m_tilde is", m_tilde)
    H_z = ustrip(H(model.cosmo, z))  # Remove units from H_z
    # println("H_z is", H_z)
    H0 = 100 * model.cosmo.h   # use h from model.cosmo
    # println("H0 is", H0)
    C = 1.65 * (model.cosmo.h / 0.7)^2 * (H_z / H0)^(8/3) *
        (m_tilde / (0.7 * 3e14))^(2/3 + 0.12)*
        (0.7 / model.cosmo.h)^(1.5)    
    # Attach the unit of eV/cm^3 to C
    C_norm = C * u"eV/cm^3"
    # println("C_norm is ", C_norm)
    return C_norm
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

#To check this A10
function dimensionless_P_profile_los(model::Arnauld10ThermalSZProfile{T}, r, M_500, z) where T
    par = get_params(model, M_500, z)
    R_500 = R_Δ(model, M_500, z, 500)/((model.B)^(1/3)) 
    x = r / angular_size(model, R_500, z)
    return par.P₀ * _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
end

"""Line-of-sight integrated electron pressure"""
P_e_los(model, r, M_200c, z) = 0.5176 * P_th_los(model, r, M_200c, z)

"""Line-of-sight integrated thermal pressure"""
P_th_los(model, r, M_200c, z) = constants.G * M_200c * 200 * ρ_crit(model, z) * 
    model.f_b / 2 * dimensionless_P_profile_los(model, r, M_200c, z)

"""
    compton_y(model, r, M_200c, z)

Calculate the Compton y parameter for a given model at a given radius, mass, and redshift.
Mass needs to have units!
"""
function compton_y(model::Battaglia16ThermalSZProfile, r, M_200c, z)
    return P_e_los(model, r, M_200c, z) * P_e_factor + 0   # +0 to strip units
end
function compton_y(model::BreakModel, r, M_200c, z)
    return P_e_los(model, r, M_200c, z) * P_e_factor + 0   # +0 to strip units
end

function compton_y(model::Arnauld10ThermalSZProfile, r, M_500, z)
    C_norm = A10_normalization(model, M_500, z; B=model.B)
    R_500 = R_Δ(model, M_500, z, 500)/((model.B)^(1/3)) 
    # println("C_norm is", C_norm)
    # println("P_e_factor is", P_e_factor)
    # println("R_200 is", R_200)
    return C_norm * dimensionless_P_profile_los(model, r, M_500, z) * R_500 * P_e_factor + 0   # +0 to strip units
end


# direct evaluation of a model at a given radius, mass, and redshift. 
# mass is in Msun; Unitful is NOT used in these evaluation functions
(model::Battaglia16ThermalSZProfile)(r, M_200c, z) = compton_y(model, r, M_200c * M_sun, z)
(model::BreakModel)(r, M_200c, z) = compton_y(model, r, M_200c * M_sun, z)
(model::Arnauld10ThermalSZProfile)(r, M_500c, z) = compton_y(model, r, M_500c * M_sun, z)
