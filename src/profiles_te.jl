"""
Implementing electron temperature (Te) moment expansion for galaxy clusters
based on SZ and tau profiles from Battaglia et al. (2016) model. 
"""

abstract type AbstractBattagliaTeMomentProfile{T,C} <: AbstractProfile{T} end

"""
    Battaglia16TeMoment0Profile{T,C} <: AbstractBattagliaTeMomentProfile{T,C}

Implements the zeroth moment of electron temperature profile, assuming Battaglia SZ and tau profiles.
Contains cosmology (C), thermal SZ profile, and tau profile.
"""
struct Battaglia16TeMoment0Profile{T,C} <: AbstractBattagliaTeMomentProfile{T,C}
    cosmo::C
    sz_model::Battaglia16ThermalSZProfile{T,C}
    tau_model::BattagliaTauProfile{T,C}
end

"""
    Battaglia16TeMoment1Profile{T,C} <: AbstractBattagliaTeMomentProfile{T,C}

Implements the first moment of electron temperature profile, assuming Battaglia SZ and tau profiles.
Contains cosmology (C), thermal SZ profile and tau profile.
"""
struct Battaglia16TeMoment1Profile{T,C} <: AbstractBattagliaTeMomentProfile{T,C}
    cosmo::C
    sz_model::Battaglia16ThermalSZProfile{T,C}
    tau_model::BattagliaTauProfile{T,C}
end

"""
    create_te_moment_profile(::Type{P}; Omega_c=0.2589, Omega_b=0.0486, h=0.6774)

Create a temperature moment profile of type P with given cosmological parameters.

# Arguments
- `P`: Profile type (Battaglia16TeMoment0Profile or Battaglia16TeMoment1Profile)
- `Omega_c`: Cold dark matter density parameter
- `Omega_b`: Baryon density parameter
- `h`: Dimensionless Hubble parameter
"""
function create_te_moment_profile(::Type{P};
    Omega_c::T=0.2589,
    Omega_b::T=0.0486,
    h::T=0.6774) where {T <: Real, P <: AbstractBattagliaTeMomentProfile}
    
    OmegaM = Omega_b + Omega_c
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    sz_model = Battaglia16ThermalSZProfile(Omega_c=Omega_c, Omega_b=Omega_b, h=h)
    tau_model = BattagliaTauProfile(Omega_c=Omega_c, Omega_b=Omega_b, h=h)
    return P(cosmo, sz_model, tau_model)
end

Battaglia16TeMoment0Profile(;kwargs...) = create_te_moment_profile(Battaglia16TeMoment0Profile; kwargs...)
Battaglia16TeMoment1Profile(;kwargs...) = create_te_moment_profile(Battaglia16TeMoment1Profile; kwargs...)

"""
    compute_meff()

Compute effective mass per electron, accounting for hydrogen and helium abundances.
Returns me + nH/ne * mH + nHe/ne * mHe.
"""
function compute_meff()
    me = constants.ElectronMass
    mH = constants.ProtonMass
    mHe = 4mH
    xH = 0.76
    nH_ne = 2xH / (xH + 1)
    nHe_ne = (1 - xH)/(2 * (1 + xH))
    return me + nH_ne*mH + nHe_ne*mHe
end

"""
    object_size(model::AbstractBattagliaTeMomentProfile, physical_size, z)

Convert physical size to angular size at redshift z using angular diameter distance.
Returns angle in radians.
"""
function object_size(model::AbstractBattagliaTeMomentProfile{T,C}, physical_size, z) where {T,C}
    d_A = angular_diameter_dist(model.cosmo, z)
    phys_siz_unitless = T(ustrip(uconvert(unit(d_A), physical_size)))
    d_A_unitless = T(ustrip(d_A))
    return atan(phys_siz_unitless, d_A_unitless)
end


"""
    T0_normalization(sz_model, tau_model, M_200c, z)

Calculate central electron temperature normalization using pressure and density profiles.
Uses SZ profile for pressure and tau profile for electron density.
"""
function T0_normalization(sz_model, tau_model, M_200c, z)
    # Compute P_e0 from SZ model
    par_sz = get_params(sz_model, M_200c, z)
    P_e0 = 0.5176 * sz_model.f_b / 2 * constants.G * M_200c * 200 * 
           par_sz.P₀ # ρ_crit(sz_model, z) cancelled

    # Compute n_e0 from τ model
    par_tau = get_params(tau_model, M_200c, z)
    R_200c = R_Δ(tau_model, M_200c, z, 200)
    m_eff = compute_meff()  # m_eff = me + nH*mH + nHe*mHe
    n_e0 = (R_200c * par_tau.P₀) / m_eff / (1+z)^2  # ρ_crit cancelled

    # T_e0 = P_e0 / (n_e0 * k_B)
    return P_e0 / (n_e0 * constants.BoltzmannConstant)
end

function dimensionless_T_profile(x, sz_model::Battaglia16ThermalSZProfile, 
    tau_model::AbstractBattagliaTauProfile, M_200c, z)
    
    par_sz = get_params(sz_model, M_200c, z)
    par_tau = get_params(tau_model, M_200c, z)
    dimensionless_T_profile(x, par_sz, par_tau)
end

function dimensionless_T_profile(x, par_sz, par_tau)
    P_gnfw = generalized_nfw(x, par_sz.xc, par_sz.α, par_sz.β, par_sz.γ)
    N_gnfw = generalized_nfw(x, par_tau.xc, par_tau.α, par_tau.β, par_tau.γ)
    return P_gnfw / N_gnfw
end

"""
    moment_los_quadrature(x, qi, sz_model, tau_model, M_200c, z; zmax=10.0, rtol=1e-4, order=7)

Compute line-of-sight integral of temperature profile raised to power (qi+1).
Uses Gauss-Kronrod quadrature for numerical integration.
"""
function moment_los_quadrature(x, qi,
    sz_model::Battaglia16ThermalSZProfile,
    tau_model::AbstractBattagliaTauProfile,
    M_200c, z::Real;
    zmax::Real=10.0, rtol::Real=1e-4, order::Integer=7)
    
    x² = x^2
    scale = 1  # Scale factor to avoid numerical underflow
    integral, err = quadgk(y -> scale * dimensionless_T_profile(√(y^2 + x²), 
        sz_model, tau_model, M_200c, z)^(qi+1),
        0.0, zmax, rtol=rtol, order=order)
    return 2integral / scale
end


"""
    Te_3d(r, model::AbstractBattagliaTeMomentProfile, M_200c, z)

Calculate 3D electron temperature profile at radius r.
Note: This function is for reference and may not be used in final implementation.
"""
function Te_3d(r, model::AbstractBattagliaTeMomentProfile, M_200c, z)
    R_200c = R_Δ(model, M_200c, z, 200)
    x = r / object_size(model, R_200c, z)  # either ang/ang or phys/phys
    T0 = T0_normalization(model.sz_model, model.tau_model, M_200c, z)
    T_dimensionless = dimensionless_T_profile(x, model.sz_model, model.tau_model, M_200c, z)
    return T0 * T_dimensionless
end

"""
    Theta_moment(model::AbstractBattagliaTeMomentProfile, qi, r, M_200c, z)

Calculate moment qi of dimensionless electron temperature θ = kB*Te/(me*c²).
Combines temperature normalization with line-of-sight integration.
"""
function Theta_moment(model::AbstractBattagliaTeMomentProfile, qi, r, M_200c, z)
    R_200c = R_Δ(model, M_200c, z, 200)
    x = r / object_size(model, R_200c, z)  # either ang/ang or phys/phys
    T_moment_dimensionless_los = moment_los_quadrature(x, qi, model.sz_model, model.tau_model, M_200c, z)
    T0 = T0_normalization(model.sz_model, model.tau_model, M_200c, z)
    factor = Theta_e_factor * T0 + 0  # ensure unitless
    return factor^(qi+1) * T_moment_dimensionless_los
end

(model::Battaglia16TeMoment0Profile)(r, M_200c, z) = 
    Theta_moment(model, 0, r, M_200c * M_sun, z)

(model::Battaglia16TeMoment1Profile)(r, M_200c, z) =
    Theta_moment(model, 1, r, M_200c * M_sun, z)