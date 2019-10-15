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

# concentration factor from Duffy et al. 2008
mz2c(m::Float64, z::Float64, cosmo) = 7.85 * (
    m / (2.0e12/cosmo.h))^(-0.081) / (1.0+z)^0.71
mz2c(m::Float32, z::Float32, cosmo) = 7.85f0 * (
    m / (2.0f12/cosmo.h))^(-0.081f0) / (1.0f0+z)^0.71f0

#
# function m2r(m, par)
#     return (3.0f0 * m / fourpif / par["rhocrit"])^(1.0f0 / 3.0f0)
# end
#
# function f_nfw(x::T) where T
#     return log1p(x) - x/(one(T)+x)
# end
# g(x,c,m) = m - f_nfw(x)/f_nfw(c)


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
