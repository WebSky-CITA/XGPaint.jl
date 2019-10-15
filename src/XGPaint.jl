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


export CIBModel




end # module
