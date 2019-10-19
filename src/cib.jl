
using Interpolations
using QuadGK
using Roots
using Cosmology
using Unitful
using UnitfulAstro
using Random
using Healpix
import Distributions
using Random

abstract type AbstractCIBModel{T<:Real} <: AbstractForegroundModel end

"""
    CIBModel_Planck2013{T}(model parameters...)

Define CIB model parameters. Defaults are from Viero et al. 2013.

```@example
model = CIBModel{Float32}(shang_Mpeak=10^12.4)
```
"""
Base.@kwdef struct CIBModel_Planck2013{T<:Real} <: AbstractCIBModel{T}
    nside::Int64    = 4096
    hod::String     = "shang"
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
    min_log_M::T, max_log_M::T, model::AbstractCIBModel;
    n_bin=1000) where T

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



function sigma_cen(m::T, model::AbstractCIBModel) where T
    return (exp( -(log10(m) - log10(model.shang_Mpeak))^2 /
        (T(2)*model.shang_sigmaM) ) * m) / sqrt(T(2π) * model.shang_sigmaM)
end

function nu2theta(nu::T, z::T, model::AbstractCIBModel) where T
    phys_h = T(6.62606957e-27)   # erg.s
    phys_k = T(1.3806488e-16)    # erg/K
    Td = model.shang_Td * (one(T)+z)^model.shang_alpha
    xnu = phys_h * nu / phys_k / Td
    return xnu^(T(4) + model.shang_beta) / expm1(xnu) / nu / model.shang_I0
end

"""
<L_sat> interpolation values
"""
function integrand_L(lm, lM_halo, model::AbstractCIBModel)
    m = exp(lm)
    return sigma_cen(m, model) * jiang_shmf(m, exp(lM_halo), model)
end


"""
Build a linear interpolator that takes in ln(M_halo) and returns sigma.
"""
function build_sigma_sat_ln_interpolator(
        max_ln_M::T, model::AbstractCIBModel; n_bin=1000) where T
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
function shang_z_evo(z::T, model::AbstractCIBModel) where T
    return (one(T) + min(z, model.shang_zplat))^model.shang_eta
end

"""
Construct the necessary interpolator set.
"""
function get_interpolators(model::AbstractCIBModel, cosmo::Cosmology.FlatLCDM{T},
    min_halo_mass::T, max_halo_mass::T) where T
    return (
        r2z = XGPaint.build_r2z_interpolator(
            model.min_redshift, model.max_redshift, cosmo),
        hod_shang = XGPaint.build_shang_interpolator(
            log(min_halo_mass), log(max_halo_mass), model),
        c_lnm2r = XGPaint.build_c_lnm2r_interpolator(),
        sigma_sat = XGPaint.build_sigma_sat_ln_interpolator(
            log(max_halo_mass), model),
        muofn = XGPaint.build_muofn_interpolator(model)
    )
end



"""
Fill up arrays with information related to CIB central sources.
"""
function process_centrals!(
    model::AbstractCIBModel{T}, cosmo::Cosmology.FlatLCDM{T}, Healpix_res::Resolution;
    interp, hp_ind_cen, dist_cen, redshift_cen,
    lum_cen, n_sat_bar, n_sat_bar_result,
    halo_pos, halo_mass) where T

    N_halos = size(halo_mass, 1)
    Threads.@threads for i = 1:N_halos
        # location information for centrals
        hp_ind_cen[i] = Healpix.vec2pixRing(Healpix_res,
            halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
        dist_cen[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
        redshift_cen[i] = interp.r2z(dist_cen[i])

        # compute HOD
        n_sat_bar[i] = interp.hod_shang(log(halo_mass[i]))
        n_sat_bar_result[i] = rand(Distributions.PoissonADSampler(
            Float64(n_sat_bar[i])))

        # get central luminosity
        lum_cen[i] = sigma_cen(halo_mass[i], model) * shang_z_evo(
            redshift_cen[i], model)
    end
end


"""
Fill up arrays with information related to CIB satellites.
"""
function process_sats!(
        model::AbstractCIBModel{T}, cosmo::Cosmology.FlatLCDM{T},
        Healpix_res::Resolution;
        interp, hp_ind_sat, dist_sat, redshift_sat,
        lum_sat, cumsat,
        halo_mass, halo_pos, redshift_cen, n_sat_bar, n_sat_bar_result) where T

    N_halos = size(halo_mass, 1)
    Threads.@threads for i_halo = 1:N_halos
        r_cen = m2r(halo_mass[i_halo], cosmo)
        c_cen = mz2c(halo_mass[i_halo], redshift_cen[i_halo], cosmo)
        for j in 1:n_sat_bar_result[i_halo]
            # i is central halo index, j is index of satellite within each halo
            i_sat = cumsat[i_halo]+j # index of satellite in satellite arrays

            log_msat_inner = max(log(rand(T)), T(-7.0))
            r_sat = r_cen * interp.c_lnm2r(
                c_cen, log_msat_inner) * T(200.0^(-1.0/3.0))
            m_sat = (interp.muofn(rand(T) * n_sat_bar[i_halo])
                * halo_mass[i_halo])

            phi = random_phi(T)
            theta = random_theta(T)
            x_sat = halo_pos[1,i_halo] + r_sat * sin(theta) * cos(phi)
            y_sat = halo_pos[2,i_halo] + r_sat * sin(theta) * sin(phi)
            z_sat = halo_pos[3,i_halo] + r_sat * cos(theta)
            dist_sat[i_sat] = sqrt(x_sat^2 + y_sat^2 + z_sat^2)
            redshift_sat[i_sat] = interp.r2z(dist_sat[i_sat])

            # lum_sat[i_sat] = interp.sigma_sat(log(m_sat)) * shang_z_evo(
            #     redshift_sat[i], model)
            lum_sat[i_sat] = sigma_cen(m_sat, model) * shang_z_evo(
                redshift_sat[i_sat], model)
            hp_ind_sat[i_sat] = Healpix.vec2pixRing(
                Healpix_res, x_sat, y_sat, z_sat)
        end
    end
end

"""
Produce a source catalog from a model and halo catalog.
"""
function generate_sources(
        # model parameters
        model::AbstractCIBModel, cosmo::Cosmology.FlatLCDM{T},
        # halo arrays
        halo_pos::Array{T,2}, halo_mass::Array{T,1};
        verbose=true) where T

    N_halos = size(halo_mass, 1)
    interp = get_interpolators( model, cosmo,
        minimum(halo_mass), maximum(halo_mass))
    res = Resolution(model.nside)

    verbose && println("Allocating for $(N_halos) centrals.")
    hp_ind_cen = Array{Int64}(undef, N_halos)  # healpix index of halo
    lum_cen = Array{T}(undef, N_halos)  # Lum of central w/o ν-dependence
    redshift_cen = Array{T}(undef, N_halos)
    dist_cen = Array{T}(undef, N_halos)
    n_sat_bar = Array{T}(undef, N_halos)
    n_sat_bar_result = Array{Int32}(undef, N_halos)

    # STEP 1: compute central properties -----------------------------------
    verbose && println("Processing centrals on $(Threads.nthreads()) threads.")
    process_centrals!(model, cosmo, res,
        interp=interp, hp_ind_cen=hp_ind_cen, dist_cen=dist_cen,
        redshift_cen=redshift_cen, lum_cen=lum_cen, n_sat_bar=n_sat_bar,
        n_sat_bar_result=n_sat_bar_result,
        halo_pos=halo_pos, halo_mass=halo_mass)

    # STEP 2: Generate satellite arrays -----------------------------
    cumsat = cumsum(n_sat_bar_result)  # set up indices for satellites
    prepend!(cumsat, 0)
    total_n_sat = cumsat[end]
    hp_ind_sat = Array{Int64}(undef, total_n_sat)  # healpix index of halo
    lum_sat = Array{T}(undef, total_n_sat)  # Lum of central w/o ν-dependence
    redshift_sat = Array{T}(undef, total_n_sat)
    dist_sat = Array{T}(undef, total_n_sat)

    # STEP 3: compute satellite properties -----------------------------------
    verbose && println("Processing $(total_n_sat) satellites.")
    process_sats!(model, cosmo, res,
        interp=interp, hp_ind_sat=hp_ind_sat, dist_sat=dist_sat,
        redshift_sat=redshift_sat,
        lum_sat=lum_sat, cumsat=cumsat,
        halo_mass=halo_mass, halo_pos=halo_pos, redshift_cen=redshift_cen,
        n_sat_bar=n_sat_bar, n_sat_bar_result=n_sat_bar_result)

    return (
        hp_ind_cen=hp_ind_cen, lum_cen=lum_cen,
        redshift_cen=redshift_cen, dist_cen=dist_cen,
        hp_ind_sat=hp_ind_sat, lum_sat=lum_sat,
        redshift_sat=redshift_sat, dist_sat=dist_sat,

        N_cen=N_halos, N_sat=total_n_sat
    )
end

"""
Paint a source catalog onto a map.
"""
function paint!(; nu_obs::T, result_map, sources, model::AbstractCIBModel) where T

    result_map .= zero(T)  # prepare the frequency map

    # process centrals for this frequency
    Threads.@threads for i in 1:sources.N_cen
        nu = (one(T) + sources.redshift_cen[i]) * nu_obs
        result_map[sources.hp_ind_cen[i]] += l2f(
            sources.lum_cen[i] * nu2theta(
                nu, sources.redshift_cen[i], model),
            sources.dist_cen[i], sources.redshift_cen[i])
    end

    # process satellites for this frequency
    Threads.@threads for i in 1:sources.N_sat
        nu = (one(T) + sources.redshift_sat[i]) * nu_obs
        result_map[sources.hp_ind_sat[i]] += l2f(
            sources.lum_sat[i] * nu2theta(
                nu, sources.redshift_sat[i], model),
            sources.dist_sat[i], sources.redshift_sat[i])
    end
end

export CIBModel_Planck2013,
    paint!,
    generate_sources,
    get_interpolators,
    build_c_lnm2r_interpolator,
    build_muofn_interpolator,
    build_r2z_interpolator,
    build_shang_interpolator,
    build_sigma_sat_ln_interpolator,
    sigma_cen
