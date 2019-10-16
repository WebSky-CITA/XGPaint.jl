
using PoissonRandom
using Interpolations
using QuadGK
using Roots
using Cosmology
using Unitful
using UnitfulAstro
using Random

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
function integrand_L(lm, lM_halo, model::CIBModel)
    m = exp(lm)
    return sigma_cen(m, model) * jiang_shmf(m, exp(lM_halo), model)
end


"""
Build a linear interpolator that takes in ln(M_halo) and returns sigma.
"""
function build_sigma_sat_ln_interpolator(
        max_ln_M::T, model; n_bin=1000) where T
    x = LinRange(log(model.shang_Mmin), max_ln_M, n_bin)
    L_mean = zero(x)

    for i in 1:n_bin
        L_mean[i], err = quadgk( t->integrand_L(t, x[i], model),
                log(model.shang_Mmin), x[i], rtol=1e-6)
    end
    return LinearInterpolation(T.(x), T.(L_mean), extrapolation_bc=zero(T))
end

function build_muofn_interpolator(model;
        min_mu::T=-6.0f0, max_mu::T=0.0f0, nbin=1000) where T
    mu = T(10) .^ LinRange(min_mu, max_mu, nbin)
    n  = zero(mu)
    for i in 1:nbin
        n[i], err = quadgk( lmu->jiang_shmf(exp(lmu), one(T), model),
                log(mu[i]), 0.0f0, rtol=1.0f-6)
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
