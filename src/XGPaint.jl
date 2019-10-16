module XGPaint

using Distributed
using HDF5
using SharedArrays
using Healpix
using PoissonRandom
using Interpolations
using QuadGK
using Base.GC
using Roots
using Cosmology
using Unitful
using UnitfulAstro
using Random

# set different seeds for worker IDs
Random.seed!(myid() + trunc(Int64, time()))


"""
    CIBModel{T}(model parameters...)

Define CIB model parameters. Defaults are from Viero et al. 2013.

```@example
model = CIBModel{Float32}(shang_Mpeak=10^12.4)
```
"""
Base.@kwdef struct CIBModel{T}
    nside::Int64    = 4096
    hod::String     = "shang"
    LM::String      = "Planck2013"
    Inu_norm::T     = 0.3180384
    min_redshift::T = 0.0
    max_redshift::T = 4.5
    min_mass::T     = 1e12
    box_size::T     = 40000

    # shang HOD
    shang_zplat::T  = 2.0
    shang_Td::T     = 20.7
    shang_beta::T   = 1.6
    shang_eta::T    = 2.4
    shang_alpha::T  = 0.2
    shang_Mpeak::T  = 10^12.3
    shang_sigmaM::T = 0.3
    shang_Msmin::T  = 1e11
    shang_Mmin::T   = 1e10
    shang_I0::T     = 46

    # jiang
    jiang_gamma_1::T    = 0.13
    jiang_alpha_1::T    = -0.83
    jiang_gamma_2::T    = 1.33
    jiang_alpha_2::T    = -0.02
    jiang_beta_2::T     = 5.67
    jiang_zeta::T       = 1.19
end

"""
Construct a cosmology with typing. The type of the output will the type of
`h` and `OmegaM`. The types of `h` and `OmegaM` must match.
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
function build_r2z_interpolator(min_z::T, max_z::T, cosmo; n_bins=1000) where T

    zrange = LinRange(min_z, max_z, n_bins)
    rrange = zero(zrange)
    for i in 1:n_bins
        rrange[i] = ustrip(T, u"Mpc",
            comoving_radial_dist(u"Mpc", cosmo, zrange[i]))
    end
    r2z = LinearInterpolation(rrange, zrange);
    return r2z
end

function jiang_shmf(m, M_halo, model)
    dndm = (((model.jiang_gamma_1*((m/M_halo)^model.jiang_alpha_1))+
             (model.jiang_gamma_2*((m/M_halo)^model.jiang_alpha_2)))*
             (exp(-(model.jiang_beta_2)*((m/M_halo)^model.jiang_zeta))))
    return dndm
end

"""
Build a linear interpolation function which maps log(M_h) to N_sat.
"""
function build_shang_interpolator(
    min_log_M::T, max_log_M::T, model::CIBModel; n_bin=1000) where T

    x_m = LinRange(min_log_M, max_log_M, 1000)
    N_sat_i = zero(x_m)

    function integrand_m(lm, lM_halo)
        # lm is natural log of m
        m = exp(lm)
        M_halo = exp(lM_halo)
        dns_dm = jiang_shmf(m, M_halo, model)
        return dns_dm
    end

    for i in 1:size(x_m,1)
        N_sat_i[i], err = quadgk( t->integrand_m(t, x_m[i]),
                log(model.shang_Msmin), x_m[i], rtol=1e-8)
        N_sat_i[i] = convert(T, max(0.0, N_sat_i[i]))
    end

    return LinearInterpolation(x_m, N_sat_i)
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

function solve_r(c,lnm)
    m = exp(lnm)
    return find_zero( x -> g(x,c,m) , (0.0,100.0), Bisection()) / c
end

"""
Generates an interpolator r(c, lnm)

Generate a LinearInterpolation object that turns concentration
and ln(M_halo) into satellite radius.
"""
function build_c_lnm2r_interp(;
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


function sigma_cen(m::T, model::CIBModel) where T
    return (exp( -(log10(m) - log10(model.shang_Mpeak))^2 /
        (T(2)*model.shang_sigmaM) ) * m) / sqrt(T(2π) * model.shang_sigmaM)
end

function nu2theta(nu::T, z::T, model::CIBModel) where T
    phys_h = T(6.62606957e-27)   # erg.s
    phys_k = T(1.3806488e-16)    # erg/K
    Td = model.shang_Td * (one(T)+z)^model.shang_alpha
    xnu = phys_h * nu / phys_k / Td
    return xnu^(T(4) + model.shang_beta) / expm1(xnu) / nu / model.shang_I0
end

"""
<L_sat> interpolation values
"""
function integrand_L(lm,lM_halo, model::CIBModel)
    m = exp(lm)
    return sigma_cen(m, model) * jiang_shmf(m, exp(lM_halo), model)
end


# need to build tests for below

"""
Build a linear interpolator that takes in ln(M_halo) and returns sigma.
"""
function build_sigma_sat_ln(min_ln_M::T, max_ln_M::T, model;
    n_bin=1000) where T
    x = LinRange(min_ln_M, max_ln_M, n_bin)
    L_mean = zero(x)

    for i in 1:n_bin
        L_mean[i], err = quadgk( t->integrand_L(t, x[i]),
                log(shang_Mmin), x[i], rtol=1e-6)
    end
    return LinearInterpolation(T.(x), T.(L_mean), extrapolation_bc=zero(T))
end

function l2f(Lum::T, r_comoving::T, redshift::T) where T
    return Lum / (T(4π) * r_comoving^2 * (one(T) + redshift) )
end

function integrand_muofn(lmu::T, model) where T
    mu = exp(lmu)
    dns_dm = jiang_shmf(mu, one(T), model)
    return dns_dm
end

function build_muofn(; min_mu::T=-6.0f0, max_mu::T=0.0f0, nbin=1000) where T
    mu = T(10) .^ LinRange(min_mu, max_mu, nbin)
    n  = zero(mu)
    for i in 1:nbin
        n[i], err = quadgk( integrand_muofn,
                log(mu[i]), 0.0f0, rtol=1.0e-6)
    end
    return LinearInterpolation(T.(reverse(n)), T.(reverse(mu)))
end

"""
Compute redshift evolution factor for LF.
"""
function shang_z_evo(z::T, model) where T
    return (one(T) + min(z, model.shang_zplat))^model.shang_eta
end


export CIBModel




# function fill_halovars!(
#         x, y, z, # inputs
#         redshift_result, dist_result) # outputs
#     """
#     This function computes distance and redshift in parallel.
#     """
#
#     N_halos = size(x,1)
#     @sync @distributed for i = 1:N_halos
#         dist_result[i] = sqrt(x[i]^2 + y[i]^2 + z[i]^2)
#         redshift_result[i] = r2z(dist_result[i])
#     end
# end
#
# function hod_shang(cen_result, sat_result, sat_bar_result,
#         halo_mass::SharedArray, par)
#     # computes shang HOD and generates a Poisson draw
#     min_lm = log(minimum(halo_mass))
#     max_lm = log(maximum(halo_mass))
#     logM_to_N = build_shang_interp(min_lm, max_lm, par)
#     N_halos = size(halo_mass,1)
#
#     @sync @distributed for i = 1:N_halos
#         cen_result[i] = 1
#
#         sat_bar_result[i] = logM_to_N(log(halo_mass[i]))
#         sat_result[i] = pois_rand(convert(Float64, sat_bar_result[i]))
#     end
#
# end


end # module
