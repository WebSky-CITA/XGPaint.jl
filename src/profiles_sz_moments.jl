"""
SZ Moment Expansion for Galaxy Clusters

This module implements moment expansion for galaxy clusters. It depends on
a pressure profile model and a electron temperature profile model. This
module specifically implements the Battaglia et al. (2016) model for electron
temperature profiles, which is derived from the tSZ and optical depth profiles.
"""

#==============================================================================#
#                              BASE TYPE                                       #
#==============================================================================#

"""
    AbstractTeProfile{T} <: AbstractProfile{T}

Abstract type for electron temperature profile models.
Concrete implementations must provide methods for calculating 
dimensionless temperature profiles.
"""
abstract type AbstractTeProfile{T} <: AbstractProfile{T} end

"""
    AbstractSZMomentProfile{T,C} <: AbstractProfile{T}

Abstract type for SZ temperature moment profile models.
These combine SZ pressure profiles with temperature profiles to 
calculate temperature moments of arbitrary order.
"""
abstract type AbstractSZMomentProfile{T,C} <: AbstractProfile{T} end

#==============================================================================#
#                           CONCRETE IMPLEMENTATIONS                           #
#==============================================================================#

"""
    Battaglia16TeProfile{T,C} <: AbstractTeProfile{T}

Electron temperature profile based on Battaglia et al. (2016) model.
Combines thermal SZ pressure profile with optical depth profile to 
derive temperature profile.

# Fields
- `cosmo::C`: Cosmological model
- `sz_model::Battaglia16ThermalSZProfile{T,C}`: Thermal SZ pressure profile
- `tau_model::BattagliaTauProfile{T,C}`: Optical depth profile

# Constructor
    Battaglia16TeProfile(; Omega_c=0.2589, Omega_b=0.0486, h=0.6774)

# Examples
```julia
# Default cosmology
te_model = Battaglia16TeProfile()

# Custom cosmology
te_model = Battaglia16TeProfile(Omega_c=0.25, Omega_b=0.05, h=0.7)
```
"""
struct Battaglia16TeProfile{T,C} <: AbstractTeProfile{T}
    cosmo::C
    sz_model::Battaglia16ThermalSZProfile{T,C}
    tau_model::BattagliaTauProfile{T,C}
end

function Battaglia16TeProfile(; 
    Omega_c::T = 0.2589, 
    Omega_b::T = 0.0486, 
    h::T = 0.6774
) where {T<:Real}
    # Validate physical parameters
    @assert 0 < Omega_c < 1 "Omega_c must be between 0 and 1, got $Omega_c"
    @assert 0 < Omega_b < 1 "Omega_b must be between 0 and 1, got $Omega_b" 
    @assert 0 < h < 2 "h must be between 0 and 2, got $h"
    @assert Omega_c + Omega_b < 1 "Total matter density Ωₘ = $(Omega_c + Omega_b) cannot exceed 1"
    
    OmegaM = Omega_b + Omega_c
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    sz_model = Battaglia16ThermalSZProfile(Omega_c=Omega_c, Omega_b=Omega_b, h=h)
    tau_model = BattagliaTauProfile(Omega_c=Omega_c, Omega_b=Omega_b, h=h)
    
    return Battaglia16TeProfile(cosmo, sz_model, tau_model)
end

"""
    Battaglia16SZMoment0Profile{T,C} <: AbstractSZMomentProfile{T,C}

Zeroth moment (qi=0) of SZ profile.
Equivalent to the standard compton y parameter profile. Simply wraps around
the Battaglia16 thermal SZ pressure profile.

# Fields
- `cosmo::C`: Cosmological model
- `sz_model::Battaglia16ThermalSZProfile{T,C}`: Thermal SZ pressure profile
"""
struct Battaglia16SZMoment0Profile{T,C} <: AbstractSZMomentProfile{T,C}
    cosmo::C
    sz_model::Battaglia16ThermalSZProfile{T,C}
end

"""
    Battaglia16SZMoment1Profile{T,C} <: AbstractSZMomentProfile{T,C}

First moment (qi=1) of temperature profile. Also known as y-weighted temperature profile.
This combines the SZ pressure profile with the electron temperature profile to
calculate the first moment of the SZ profile.

# Fields
- `cosmo::C`: Cosmological model
- `sz_model::Battaglia16ThermalSZProfile{T,C}`: Thermal SZ pressure profile
- `te_model::Battaglia16TeProfile{T,C}`: Temperature profile model
"""
struct Battaglia16SZMoment1Profile{T,C} <: AbstractSZMomentProfile{T,C}
    cosmo::C
    sz_model::Battaglia16ThermalSZProfile{T,C}
    te_model::Battaglia16TeProfile{T,C}
end

"""
    Battaglia16SZMoment2Profile{T,C} <: AbstractSZMomentProfile{T,C}

Second moment (qi=2) of SZ temperature profile.

# Fields
- `cosmo::C`: Cosmological model
- `sz_model::Battaglia16ThermalSZProfile{T,C}`: Thermal SZ pressure profile
- `te_model::Battaglia16TeProfile{T,C}`: Temperature profile model
"""
struct Battaglia16SZMoment2Profile{T,C} <: AbstractSZMomentProfile{T,C}
    cosmo::C
    sz_model::Battaglia16ThermalSZProfile{T,C}
    te_model::Battaglia16TeProfile{T,C}
end

"""
    create_sz_moment_profile(qi::Integer; kwargs...) -> AbstractSZMomentProfile

Create a temperature moment profile for moment order qi.

# Arguments
- `qi`: Moment order (0, 1, or 2)
- `Omega_c`: Cold dark matter density parameter (default: 0.2589)
- `Omega_b`: Baryon density parameter (default: 0.0486)
- `h`: Dimensionless Hubble parameter (default: 0.6774)

# Returns
- `Battaglia16SZMoment0Profile` for qi = 0
- `Battaglia16SZMoment1Profile` for qi = 1  
- `Battaglia16SZMoment2Profile` for qi = 2

# Examples
```julia
# Create different moment profiles
moment0 = create_sz_moment_profile(0)
moment1 = create_sz_moment_profile(1, h=0.7)
moment2 = create_sz_moment_profile(2, Omega_c=0.25)
```
"""
function create_sz_moment_profile(
    qi::Integer;
    Omega_c::T = 0.2589,
    Omega_b::T = 0.0486,
    h::T = 0.6774
) where {T<:Real}
    
    if !(qi in [0, 1, 2])
        throw(ArgumentError("Moment order qi must be 0, 1, or 2, got $qi"))
    end
    
    OmegaM = Omega_b + Omega_c
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    sz_model = Battaglia16ThermalSZProfile(Omega_c=Omega_c, Omega_b=Omega_b, h=h)
    
    if qi == 0
        return Battaglia16SZMoment0Profile(cosmo, sz_model)
    else  # qi == 1 or 2
        te_model = Battaglia16TeProfile(Omega_c=Omega_c, Omega_b=Omega_b, h=h)
        if qi == 1
            return Battaglia16SZMoment1Profile(cosmo, sz_model, te_model)
        else  # qi == 2
            return Battaglia16SZMoment2Profile(cosmo, sz_model, te_model)
        end
    end
end

# Convenience constructors
Battaglia16SZMoment0Profile(; kwargs...) = create_sz_moment_profile(0; kwargs...)
Battaglia16SZMoment1Profile(; kwargs...) = create_sz_moment_profile(1; kwargs...)
Battaglia16SZMoment2Profile(; kwargs...) = create_sz_moment_profile(2; kwargs...)


#==============================================================================#
#                          UTILITIES                                           #
#==============================================================================#

"""
    compute_meff(; xH::Real = 0.76) -> Real

Compute effective mass per electron accounting for hydrogen and helium abundances.

The effective mass includes contributions from electrons, protons, and helium nuclei:
    mₑff = mₑ + (nH/nₑ) × mH + (nHe/nₑ) × mHe

# Arguments
- `xH`: Hydrogen mass fraction (default: 0.76, primordial value)

# Returns
- Effective mass per electron in appropriate units

# Physics
For a fully ionized plasma with hydrogen fraction xH:
- nH/nₑ = 2xH/(xH + 1)
- nHe/nₑ = (1 - xH)/(2(1 + xH))
"""
function compute_meff(; xH::Real = 0.76)
    me = constants.ElectronMass
    mH = constants.ProtonMass
    mHe = 4 * mH
    
    # Number density ratios for fully ionized plasma
    nH_ne = 2 * xH / (xH + 1)
    nHe_ne = (1 - xH) / (2 * (1 + xH))
    
    return me + nH_ne * mH + nHe_ne * mHe
end

"""
    object_size(model::AbstractSZMomentProfile, physical_size, z) -> Real

Convert physical size to angular size at redshift z.

Uses the angular diameter distance to convert from physical to angular coordinates.

# Arguments
- `model`: SZ moment profile model (contains cosmology)
- `physical_size`: Physical size (with units)
- `z`: Redshift

# Returns
- Angular size in radians

# Physics
Angular size θ = D_physical / D_A(z) where D_A is the angular diameter distance.
"""
function object_size(
    model::AbstractSZMomentProfile{T,C}, 
    physical_size, 
    z::Real
) where {T,C}
    d_A = angular_diameter_dist(model.cosmo, z)
    phys_size_unitless = T(ustrip(uconvert(unit(d_A), physical_size)))
    d_A_unitless = T(ustrip(d_A))
    return atan(phys_size_unitless, d_A_unitless)
end

#==============================================================================#
#                           PROFILE CALCULATIONS                              #
#==============================================================================#

"""
    dimensionless_T_profile(model::Battaglia16TeProfile, x::Real, M_200c::Real, z::Real) -> Real

Calculate dimensionless temperature profile T̃(x) = P(x)/N(x).

This combines the generalized NFW pressure profile P(x) with the 
number density profile N(x) to obtain the temperature structure.

# Arguments
- `model`: Temperature profile model
- `x`: Dimensionless radius (r/R₂₀₀c)
- `M_200c`: Cluster mass at 200× critical density [M☉]
- `z`: Redshift

# Returns
- Dimensionless temperature ratio T̃(x)

# Physics
The temperature profile is derived from the ideal gas law:
    T(x) = P(x)/(nₑ(x) × kB) ∝ P(x)/N(x)
where P(x) and N(x) are generalized NFW profiles with different parameters.
"""
function dimensionless_T_profile(
    model::Battaglia16TeProfile,
    x, M_200c, z
)
    par_sz = get_params(model.sz_model, M_200c, z)
    par_tau = get_params(model.tau_model, M_200c, z)
    
    P_gnfw = generalized_nfw(x, par_sz.xc, par_sz.α, par_sz.β, par_sz.γ)
    N_gnfw = generalized_nfw(x, par_tau.xc, par_tau.α, par_tau.β, par_tau.γ)
    
    return P_gnfw / N_gnfw
end

"""
    Theta_normalization(model::Battaglia16TeProfile, M_200c::Real, z::Real) -> Real

Computes the dimensionless temperature Θ₀ = kBTₑ,₀/(mₑc²) at the cluster center
by combining pressure and density normalizations from the SZ and τ profiles.

# Arguments
- `model`: Battaglia16 temperature profile model
- `M_200c`: Cluster mass at 200× critical density [M☉]
- `z`: Redshift

"""
function Theta_normalization(
    model::Battaglia16TeProfile, 
    M_200c, z
)
    # Compute electron pressure normalization from SZ model
    par_sz = get_params(model.sz_model, M_200c, z)
    P_e0 = (0.5176 * model.sz_model.f_b / 2 * constants.G * M_200c * 200 * 
            par_sz.P₀) # ρ_crit(sz_model, z) terms cancel

    # Compute electron number density normalization from τ model
    par_tau = get_params(model.tau_model, M_200c, z)
    R_200c = R_Δ(model.tau_model, M_200c, z, 200)
    m_eff = compute_meff()
    n_e0 = (R_200c * par_tau.P₀) / (m_eff * (1 + z)^2) # ρ_crit terms cancel
    
    # Central temperature
    T_e0 = P_e0 / (n_e0 * constants.BoltzmannConstant)
    
    # Convert to dimensionless Theta
    return Theta_e_factor * T_e0 + 0  # ensure unitless
end

"""
    Theta_3d(model::AbstractTeProfile, x::Real, M_200c::Real, z::Real) -> Real

Calculate 3D electron temperature profile at dimensionless radius x.

Combines the central temperature normalization with the dimensionless 
temperature profile to give the full 3D temperature structure.

# Arguments
- `model`: Temperature profile model
- `x`: Dimensionless radius (r/R₂₀₀c)
- `M_200c`: Cluster mass at 200× critical density [M☉]
- `z`: Redshift

# Returns
- Dimensionless temperature Θ(x) = kBTₑ(x)/(mₑc²)

# Examples
```julia
model = Battaglia16TeProfile()
theta_center = Theta_3d(model, 0.0, 1e14, 0.3)  # Center
theta_r200 = Theta_3d(model, 1.0, 1e14, 0.3)    # At R₂₀₀c
```
"""
function Theta_3d(
    model::Battaglia16TeProfile,
    x,
    M_200c,
    z
)
    Theta0 = Theta_normalization(model, M_200c, z)
    T_dimensionless = dimensionless_T_profile(model, x, M_200c, z)
    return Theta0 * T_dimensionless
end

# LOS integration
"""
    sz_moment_los_quadrature(x, qi, model, M_200c, z; kwargs...) -> Real

Compute line-of-sight integral of temperature profile raised to power qi.

Performs the integral:
    ∫₀^∞ P(√(x² + y²)) × Θ(√(x² + y²))^qi dy

Uses adaptive Gauss-Kronrod quadrature for numerical integration.

# Arguments
- `x`: Dimensionless projected radius
- `qi`: Temperature moment order (0, 1, or 2)
- `model`: SZ moment profile model
- `M_200c`: Cluster mass at 200× critical density [M☉]
- `z`: Redshift

# Keyword Arguments
- `zmax`: Maximum integration distance (default: 10.0)
- `rtol`: Relative tolerance for integration (default: 1e-4)
- `order`: Quadrature order (default: 7)

# Returns
- Line-of-sight integrated moment value

# Notes
- Uses scaling factor of 1e9 to avoid numerical underflow
- For qi > 0, requires temperature model in the profile struct
- Integration performed in cylindrical coordinates
"""
function sz_moment_los_quadrature(
    x,
    qi::Integer,
    model::AbstractSZMomentProfile,
    M_200c,
    z,
    zmax::Real = 10.0, 
    rtol::Real = 1e-4, 
    order::Integer = 7
)
    @assert qi in [0, 1, 2] "Moment order qi must be 0, 1, or 2, got $qi"

    # Get SZ profile parameters
    par_sz = get_params(model.sz_model, M_200c, z)
    x² = x^2
    scale = 1e9  # Scale factor to avoid numerical underflow
    
    # Check for temperature model requirement
    if qi > 0 && !hasfield(typeof(model), :te_model)
        throw(ArgumentError("Temperature model required for qi > 0"))
    end
    
    if qi > 0
        integral, err = quadgk(
            y -> scale * generalized_nfw(√(y^2 + x²), par_sz.xc, par_sz.α, par_sz.β, par_sz.γ) *
                 Theta_3d(model.te_model, √(y^2 + x²), M_200c, z)^qi,
            0.0, zmax, rtol=rtol, order=order)
    else
        @show x
        integral, err = quadgk(
            y -> scale * generalized_nfw(√(y^2 + x²), par_sz.xc, par_sz.α, par_sz.β, par_sz.γ),
            0.0, zmax, rtol=rtol, order=order)
    end
    return 2integral / scale
end

#==============================================================================#
#                          MAIN INTERFACE FOR SZ MOMENTS                       #
#==============================================================================#

"""
    sz_moment(model::AbstractSZMomentProfile, qi::Integer, r, M_200c, z) -> Real

Calculate moment qi of dimensionless electron temperature θ = kBTₑ/(mₑc²).

This is the main function that combines temperature normalization with 
line-of-sight integration to compute observable SZ temperature moments.

# Arguments
- `model`: SZ moment profile model
- `qi`: Moment order (0, 1, or 2)
- `r`: Radius (angular, in radians)
- `M_200c`: Cluster mass at 200× critical density [M☉]
- `z`: Redshift

# Returns
- Temperature moment ⟨Θⁱ⟩(r) [dimensionless]

# Physics
Calculates the projected moment:
    ⟨Θⁱ⟩(R) = ∫ Θ(r)ⁱ × nₑ(r) × dl
where the integration is along the line of sight at projected radius R.

# Examples
```julia
model1 = Battaglia16SZMoment1Profile()
first_moment = sz_moment(model1, 1, 0.001, 1e14, 0.3)  # First moment at small radius
```
"""
function sz_moment(
    model::AbstractSZMomentProfile, 
    qi::Integer, 
    r, M_200c, z
)
    sz_model = model.sz_model
    par_sz = get_params(sz_model, M_200c, z)
    R_200c = R_Δ(sz_model, M_200c, z, 200)
    
    # Convert radius to dimensionless form
    x = r / object_size(model, R_200c, z)
    
    # Calculate dimensionless line-of-sight moment
    sz_moment_dimensionless_los = sz_moment_los_quadrature(x, qi, model, M_200c, z)
    
    # Apply physical prefactor
    factor = (P_e_factor * 0.5176 * constants.G * M_200c * 200 * 
              ρ_crit(sz_model, z) * sz_model.f_b / 2 * par_sz.P₀)
    
    return factor * sz_moment_dimensionless_los + 0  # ensure unitless
end

# Callable interface for profile structs
(model::Battaglia16SZMoment0Profile)(r, M_200c, z) = 
    sz_moment(model, 0, r, M_200c * M_sun, z)

(model::Battaglia16SZMoment1Profile)(r, M_200c, z) =
    sz_moment(model, 1, r, M_200c * M_sun, z)

(model::Battaglia16SZMoment2Profile)(r, M_200c, z) =
    sz_moment(model, 2, r, M_200c * M_sun, z)