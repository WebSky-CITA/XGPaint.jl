"""Basic functions for computing models."""

"""
All foreground models inherit from this type.
"""
abstract type AbstractForegroundModel end

"""
Construct a background cosmology.

This function duplicates the cosmology() function in Cosmology.jl, but with
typing. The type of the cosmology will the type of `h` and `OmegaM`. This is
primarily for keeping the code entirely in Float32 or Float64.
"""
function get_cosmology(;h::T=0.69,
                   Neff=3.04,
                   OmegaK=0.0,
                   OmegaM::T=0.29,
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


"""
Construct a fast r2z linear interpolator.
"""
function build_r2z_interpolator(min_z::T, max_z::T,
    cosmo::Cosmology.AbstractCosmology; n_bins=1000) where T

    zrange = LinRange(min_z, max_z, n_bins)
    rrange = zero(zrange)
    for i in 1:n_bins
        rrange[i] = ustrip(T, u"Mpc",
            Cosmology.comoving_radial_dist(u"Mpc", cosmo, zrange[i]))
    end
    r2z = LinearInterpolation(rrange, zrange);
    return r2z
end

"""
Compute concentration factor from Duffy et al. 2008.
"""
function mz2c(m::T, z::T, cosmo::Cosmology.FlatLCDM{T}) where T
    return T(7.85) * (m / (T(2.0e12)/cosmo.h))^(T(-0.081)) / (T(1.0)+z)^T(0.71)
end

"""
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
Generates an interpolator r(c, lnm)

Generate a LinearInterpolation object that turns concentration
and ln(M_halo) into satellite radius.
"""
function build_c_lnm2r_interpolator(;
    cmin::T=1f-3, cmax::T=25.0f0, mmin::T=-7.1f0, mmax::T=0.0f0,
    nbin=100) where T
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


export get_cosmology, m2r, mz2c, random_phi, random_theta, l2f
