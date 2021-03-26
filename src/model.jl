"""Basic functions for computing models."""

"""
All foreground models inherit from this type.
"""
abstract type AbstractForegroundModel end

"""
    get_cosmology(::Type{T}; h=0.69, Neff=3.04, OmegaK=0.0,
        OmegaM=0.29, OmegaR=nothing, Tcmb=2.7255, w0=-1, wa=0)

Construct a background cosmology. This function duplicates the cosmology() function
in Cosmology.jl, but with typing. This is primarily for keeping the code entirely
in Float32 or Float64.

# Arguments:
- `::Type{T}`: numerical type to use for calculations

# Keywords
- `h` - Dimensionless Hubble constant
- `OmegaK` - Curvature density (Ω_k)
- `OmegaM` - Matter density (Ω_m)
- `OmegaR` - Radiation density (Ω_r)
- `Tcmb` - CMB temperature in Kelvin; used to compute Ω_γ
- `Neff` - Effective number of massless neutrino species; used to compute Ω_ν
- `w0` - CPL dark energy equation of state; `w = w0 + wa(1-a)`
- `wa` - CPL dark energy equation of state; `w = w0 + wa(1-a)`

# Example
```julia-repl
julia> get_cosmology(Float32; h=0.7)
Cosmology.FlatLCDM{Float32}(0.7f0, 0.7099147f0, 0.29f0, 8.5307016f-5)
```
"""
function get_cosmology(::Type{T}; h=0.69,
                   Neff=3.04,
                   OmegaK=0.0,
                   OmegaM=0.29,
                   OmegaR=nothing,
                   Tcmb=2.7255,
                   w0=-1,
                   wa=0) where T

    if OmegaR === nothing
        OmegaG = 4.48131e-7*Tcmb^4/h^2
        OmegaN = Neff*OmegaG*(7/8)*(4/11)^(4/3)
        OmegaR = OmegaG + OmegaN
    end

    OmegaL = 1 - OmegaK - OmegaM - OmegaR

    if !(w0 == -1 && wa == 0)
        return Cosmology.WCDM{T}(h, OmegaK, OmegaL, OmegaM, OmegaR, w0, wa)
    end

    if OmegaK < 0
        return Cosmology.ClosedLCDM{T}(h, OmegaK, OmegaL, OmegaM, OmegaR)
    elseif OmegaK > 0
        return Cosmology.OpenLCDM{T}(h, OmegaK, OmegaL, OmegaM, OmegaR)
    else
        return Cosmology.FlatLCDM{T}(h, OmegaL, OmegaM, OmegaR)
    end
end
get_cosmology(; h=0.69, Neff=3.04, OmegaK=0.0, OmegaM=0.29, OmegaR=nothing, Tcmb=2.7255,
    w0=-1, wa=0) where T = get_cosmology(Float32; h=h, Neff=Neff, OmegaK=OmegaK,
        OmegaM=OmegaM, OmegaR=OmegaR, Tcmb=Tcmb, w0=w0, wa=wa)

"""
    build_r2z_interpolator(min_z::T, max_z::T,
        cosmo::Cosmology.AbstractCosmology; n_bins=2000) where T

Construct a fast r2z linear interpolator.
"""
function build_r2z_interpolator(min_z::T, max_z::T,
    cosmo::Cosmology.AbstractCosmology; n_bins=2000) where T

    zrange = LinRange(min_z, max_z, n_bins)
    rrange = zero(zrange)
    for i in 1:n_bins
        rrange[i] = ustrip(T, u"Mpc",
            Cosmology.comoving_radial_dist(u"Mpc", cosmo, zrange[i]))
    end
    r2z = LinearInterpolation(rrange, zrange; extrapolation_bc=Line());
    return r2z
end

"""
    mz2c(m::T, z::T, cosmo::Cosmology.FlatLCDM{T}) where T

Compute concentration factor from Duffy et al. 2008.
"""
function mz2c(m::T, z::T, cosmo::Cosmology.FlatLCDM{T}) where T
    return T(7.85) * (m / (T(2.0e12)/cosmo.h))^(T(-0.081)) / (T(1.0)+z)^T(0.71)
end

"""
    m2r(m::T, cosmo::Cosmology.FlatLCDM{T}) where T

Convert virial mass to virial radius.
"""
function m2r(m::T, cosmo::Cosmology.FlatLCDM{T}) where T
    rho = T(2.78e11) * cosmo.Ω_m * cosmo.h^2
    return (T(3.0/(4π)) * m / rho)^T(1.0/3.0)
end


function f_nfw(x::T) where T
    return log1p(x) - x/(one(T)+x)
end
g(x, c, m) = m - f_nfw(x)/f_nfw(c)

function solve_r(c::T, lnm::T) where T
    m = exp(lnm)
    return find_zero( x -> g(x,c,m) , (T(0.0),T(100.0)), Bisection()) / c
end


"""
    build_c_lnm2r_interpolator(;cmin::T=1f-3, cmax::T=25.0f0,
        mmin::T=-7.1f0, mmax::T=0.0f0, nbin=100) where T

Generate a LinearInterpolation object that turns concentration
and ln(M_halo) into satellite radius.
"""
function build_c_lnm2r_interpolator(;cmin::T=1f-3, cmax::T=25.0f0,
        mmin::T=-7.1f0, mmax::T=0.0f0, nbin=100) where T
    # these defaults imply a minimum fractional mass of exp(-7)

    cs  = LinRange(cmin,cmax,nbin)
    ms  = LinRange(mmin,mmax,nbin)
    r = [solve_r(c_,m_) for c_ in cs, m_ in ms]

    return LinearInterpolation((T.(cs), T.(ms)), T.(r))
end

random_phi(::Type{Float32}) = rand(Float32) * Float32(2π)
random_theta(::Type{Float32}) = acos(2.0f0*rand(Float32)-1.0f0)
random_phi(t::DataType) = rand(t) * t(2π)
random_theta(t::DataType) = acos( t(2.0)*rand(t)-one(t))

"""
Inverse square law with redshift dependence.
"""
function l2f(luminosity::T, r_comoving::T, redshift::T) where T
    return luminosity / (T(4π) * r_comoving^2 * (one(T) + redshift) )
end


function flux2map!(result_map::Map{T_map,RingOrder}, fluxes, theta, phi) where {T_map, T}

    pixel_array = result_map.pixels
    fill!(pixel_array, zero(T))  # prepare the frequency map

    flux_to_Jy = ustrip(u"Jy", 1u"W/Hz*Mpc^-2")

    source_offset_I = generate_subhalo_offsets(sources.nsources_I)
    source_offset_II = generate_subhalo_offsets(sources.nsources_II)
    N_halo = size(sources.halo_mass, 1)
    total_n_sat_I = size(sources.L_I_151, 1)
    total_n_sat_II = size(sources.L_II_151, 1)

    flux_I = Array{T, 1}(undef, total_n_sat_I)
    flux_II = Array{T, 1}(undef, total_n_sat_II)
    redshift_I = Array{T, 1}(undef, total_n_sat_I)
    redshift_II = Array{T, 1}(undef, total_n_sat_II)
    θ_I = Array{T, 1}(undef, total_n_sat_I)
    θ_II = Array{T, 1}(undef, total_n_sat_II)
    ϕ_I = Array{T, 1}(undef, total_n_sat_I)
    ϕ_II = Array{T, 1}(undef, total_n_sat_II)


    Threads.@threads for i_halo in 1:N_halo
        nu = (1 + sources.redshift[i_halo]) * nu_obs
        hp_ind = sources.halo_hp_ind[i_halo]

        # Generate FR I fluxes
        for source_index in 1:sources.nsources_I[i_halo]
            # index of satellite in satellite arrays
            i_sat = source_offset_I[i_halo] + source_index
            L_c, L_l = get_core_lobe_lum(
                sources.L_I_151[i_sat], nu, model.I_R_int, model.I_γ,
                model.a_0,
                sources.a_coeff_I[1,i_sat], sources.a_coeff_I[2,i_sat],
                sources.a_coeff_I[3,i_sat], sources.cosθ_I[i_sat])

            flux_I[i_sat] = l2f(L_c+L_l, sources.dist[i_halo],
                sources.redshift[i_halo]) * flux_to_Jy
            redshift_I[i_sat] = sources.redshift[i_halo]
            θ_I[i_sat] = sources.θ[i_halo]
            ϕ_I[i_sat] = sources.ϕ[i_halo]
            pixel_array[hp_ind] += flux_I[i_sat]
        end

        # Generate FR II fluxes
        for source_index in 1:sources.nsources_II[i_halo]
            # index of satellite in satellite arrays
            i_sat = source_offset_II[i_halo] + source_index
            L_c, L_l = get_core_lobe_lum(
                sources.L_II_151[i_sat], nu, model.II_R_int, model.II_γ,
                model.a_0,
                sources.a_coeff_II[1,i_sat], sources.a_coeff_II[2,i_sat],
                sources.a_coeff_II[3,i_sat], sources.cosθ_II[i_sat])

            flux_II[i_sat] = l2f(L_c+L_l, sources.dist[i_halo],
                sources.redshift[i_halo]) * flux_to_Jy
            redshift_II[i_sat] = sources.redshift[i_halo]
            θ_II[i_sat] = sources.θ[i_halo]
            ϕ_II[i_sat] = sources.ϕ[i_halo]
            pixel_array[hp_ind] += flux_II[i_sat]
        end
    end

    # divide by healpix pixel size
    per_pixel_steradian = 1 / nside2pixarea(result_map.resolution.nside)
    pixel_array .*= per_pixel_steradian

    # maps are in Jansky per steradian, fluxes are in Jansky
    if return_fluxes
        return flux_I, redshift_I, θ_I, ϕ_I, flux_II, redshift_II, θ_II, ϕ_II
    end

end
