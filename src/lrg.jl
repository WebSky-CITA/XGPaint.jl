"""Functions for computing a very rudimentary DESI LRG mock."""

abstract type AbstractLRGModel{T<:Real} <: AbstractForegroundModel end

"""
    LRG_Yuan22{T}(; kwargs...)

Define LRG model parameters. Defaults from Yuan et al. 2022 [arXiv:2202.12911]. All numbers
not typed are converted to type T. This model has the following parameters and default values:

* `nside::Int64 = 4096`
* `hod::String = "zheng"`
* `min_redshift = 0.0`
* `max_redshift = 5.0`
* `min_mass = 1e12`
* `box_size = 40000`
* `shang_Msmin = 1e11`
* `zheng_Mcut = 10^12.7`
* `zheng_M1 = 10^13.6`
* `zheng_sigma = 0.2`
* `zheng_kappa = 0.08`
* `zheng_alpha = 1.15`
* `yuan_ic = 0.8`
* `jiang_gamma_1 = 0.13`
* `jiang_alpha_1 = -0.83`
* `jiang_gamma_2 = 1.33`
* `jiang_alpha_2 = -0.02`
* `jiang_beta_2 = 5.67`
* `jiang_zeta = 1.19`

"""
@with_kw struct LRG_Yuan22{T<:Real} <: AbstractLRGModel{T} @deftype T
    nside::Int64    = 4096
    hod::String     = "zheng"
    min_redshift = 0.0
    max_redshift = 5.0
    min_mass     = 1e12
    box_size     = 40000
    
    # need to account for where this param is actually from
    shang_Msmin = 1e11
    
    # UniverseMachine-derived quenching fraction
    quench_account::Bool = false
    quench_Qmin0 = -1.944
    quench_Qmina = -2.419
    quench_VQ0 = 2.248
    quench_VQa = 0.018
    quench_VQz = 0.124
    quench_sigVQ0 = 0.227
    quench_sigVQa = 0.037
    quench_sigVQl = 0.107
    fquench_max  = 1.0

    # zheng HOD
    zheng_Mcut = 10^12.7/0.677
    zheng_M1 = 10^13.6/0.677
    zheng_sigma = 0.2
    zheng_kappa = 0.08
    zheng_alpha = 1.15
    yuan_ic = 0.8

    # jiang
    jiang_gamma_1    = 0.13
    jiang_alpha_1    = -0.83
    jiang_gamma_2    = 1.33
    jiang_alpha_2    = -0.02
    jiang_beta_2     = 5.67
    jiang_zeta       = 1.19
end

"""
[already in cib.jl, reproduced for reference]
function jiang_shmf(m, M_halo, model)
    dndm = (((model.jiang_gamma_1*((m/M_halo)^model.jiang_alpha_1))+
             (model.jiang_gamma_2*((m/M_halo)^model.jiang_alpha_2)))*
             (exp(-(model.jiang_beta_2)*((m/M_halo)^model.jiang_zeta))))
    return dndm
end
"""

"""
Build a linear interpolation function which maps log(M_h) to N_cen.
"""
function build_zhengcen_interpolator(
    min_log_M::T, max_log_M::T, model::AbstractLRGModel;
    n_bin=1000) where T

    x_m = LinRange(min_log_M, max_log_M, 1000)
    N_cen_i = zero(x_m)
    l10Mcut = log10(model.zheng_Mcut)

    for i in 1:size(x_m,1)
        N_cen_i[i] = model.yuan_ic/T(2)*SpecialFunctions.erfc((l10Mcut-x_m[i]/2.30258509)/(T(1.41421356)*model.zheng_sigma))
        N_cen_i[i] = convert(T, max(0.0, N_cen_i[i]))
    end

    return LinearInterpolation(x_m, N_cen_i)
end

"""
Build a linear interpolation function which maps log(M_h) to N_sat.
"""
function build_zhengsat_interpolator(
    min_log_M::T, max_log_M::T, model::AbstractLRGModel;
    n_bin=1000) where T

    x_m = LinRange(min_log_M, max_log_M, 1000)
    N_sat_i = zero(x_m)

    for i in 1:size(x_m,1)
        N_sat_i[i] = max(T(0),(exp(x_m[i])-model.zheng_kappa*model.zheng_Mcut)/model.zheng_M1)^model.zheng_alpha
        N_sat_i[i] = convert(T, max(0.0, N_sat_i[i]))
    end

    return LinearInterpolation(x_m, N_sat_i)
end

"""
Build a linear interpolation function which maps log(M_h) to N_sh.
"""
function build_jiang_interpolator(
    min_log_M::T, max_log_M::T, model::AbstractLRGModel;
    n_bin=1000) where T

    x_m = LinRange(min_log_M, max_log_M, 1000)
    N_sh_i = zero(x_m)

    function integrand_m(lm, lM_halo)
        # lm is natural log of m
        m = exp(lm)
        M_halo = exp(lM_halo)
        dns_dm = jiang_shmf(m, M_halo, model)
        return dns_dm
    end

    for i in 1:size(x_m,1)
        N_sh_i[i], err = quadgk( t->integrand_m(t, x_m[i]),
                log(model.shang_Msmin), x_m[i], rtol=1e-8)
        N_sh_i[i] = convert(T, max(0.0, N_sh_i[i]))
    end

    return LinearInterpolation(x_m, N_sh_i)
end

"""
Build a linear interpolator for drawing masses from the subhalo mass function.
[already in cib.jl, avoiding duplication]

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

"""
Quiescent fraction recipe from UniverseMachine
"""
function fquench_UM(Mh::T,z::T,model::AbstractLRGModel) where T
    a = one(T)/(one(T)+z);
    M200kms = T(1.64e12)/((a/T(0.378))^T(-0.142)+(a/T(0.378))^T(-1.79)) # MSol
    vMpeak = T(200)*(Mh/M200kms)^T(0.3) # km/s
    Qmin = max(T(0),T(model.quench_Qmin0)+T(model.quench_Qmina)*(one(T)-a));
    logVQ = T(model.quench_VQ0) - T(model.quench_VQa)*(T(1)-a) + T(model.quench_VQz)*z;
    sigVQ = T(model.quench_sigVQ0) + T(model.quench_sigVQa)*(T(1)-a) - T(model.quench_sigVQl)*log(T(1)+z);
    return Qmin + (one(T)-Qmin)*(T(0.5)+T(0.5)*erf((log10(vMpeak)-logVQ)/(T(sqrt(2))*sigVQ)))
end

function build_fquench_interpolator(
    max_log_M::T, model::AbstractLRGModel;
    n_bin=1000) where T

    logMh_range = LinRange(log(model.shang_Msmin), max_log_M, n_bin)
    log1z_range = LinRange(log(one(T)+model.min_redshift),log(one(T)+model.max_redshift),n_bin)

    fquench_table = [fquench_UM(exp(Mh),exp(zp1)-1,model) for Mh in logMh_range, zp1 in log1z_range]

    return linear_interpolation((logMh_range,log1z_range),fquench_table)
end


"""
Construct the necessary interpolator set.
"""
function get_interpolators(model::AbstractLRGModel, cosmo::Cosmology.FlatLCDM{T},
    min_halo_mass::T, max_halo_mass::T) where T
    return (
        r2z = build_r2z_interpolator(
            model.min_redshift, model.max_redshift, cosmo),
        hod_jiang = build_jiang_interpolator(
            log(min_halo_mass), log(max_halo_mass), model),
        hod_zhengcen = build_zhengcen_interpolator(
            log(min_halo_mass), log(max_halo_mass), model),
        hod_zhengsat = build_zhengsat_interpolator(
            log(min_halo_mass), log(max_halo_mass), model),
        c_lnm2r = build_c_lnm2r_interpolator(),
        muofn = build_muofn_interpolator(model),
        fquench = build_fquench_interpolator(
            log(max_halo_mass), model)
    )
end



# Fill up arrays with information related to CIB central sources.
function process_centrals!(
    model::AbstractLRGModel{T}, cosmo::Cosmology.FlatLCDM{T}, Healpix_res::Resolution;
    interp, hp_ind_cen, dist_cen, redshift_cen, theta_cen, phi_cen, m200c_cen,
    lrg_cen, n_sh_bar, n_sh_bar_result,
    halo_pos, halo_mass, halo_quenched) where T

    N_halos = size(halo_mass, 1)

    Threads.@threads for i = 1:N_halos
        # location information for centrals
        hp_ind_cen[i] = Healpix.vec2pixRing(Healpix_res,
            halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
        theta_cen[i], phi_cen[i] = Healpix.vec2ang(halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
        dist_cen[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
        redshift_cen[i] = interp.r2z(dist_cen[i])
        # deformation of Bocquet+16
        m200c_cen[i] = halo_mass[i]*((0.01*redshift_cen[i]-0.0225)*log10(halo_mass[i])-0.0197*redshift_cen[i]+1.0564)
        # compute central HOD
        lrg_frac = interp.hod_zhengcen(log(m200c_cen[i]))
        if model.quench_account
            lrg_cen[i] = false
            fquench_result = min(model.fquench_max,interp.fquench(log(halo_mass[i]),log(1+redshift_cen[i])))
            if fquench_result >= lrg_frac
                if halo_quenched[i]
                    lrg_cen[i] = rand(T) < lrg_frac/fquench_result
                end
            else
                if halo_quenched[i]
                    lrg_cen[i] = true
                else
                    lrg_cen[i] = rand(T) < (lrg_frac-fquench_result)/(1-fquench_result)
                end
            end
        else
            lrg_cen[i] = rand(T) < lrg_frac
        end
        # compute ***total*** subhalo (not satellite LRG) count
        n_sh_bar[i] = interp.hod_jiang(log(m200c_cen[i]))
        n_sh_bar_result[i] = rand(Distributions.Poisson(Float64.(n_sh_bar[i])))
    end
end


# Fill up arrays with information related to CIB satellites.
function process_sats!(
        model::AbstractLRGModel{T}, cosmo::Cosmology.FlatLCDM{T},
        Healpix_res::Resolution;
        interp, hp_ind_sh, dist_sh, parent_sh, redshift_sh, theta_sh, phi_sh,
        m200c_sh, lrg_sh, cumush, lrg_cen,
        m200c_cen, halo_pos, redshift_cen, n_sh_bar, n_sh_bar_result) where T

    N_halos = size(m200c_cen, 1)
    Threads.@threads for i_halo = 1:N_halos
        r_cen = m2r(m200c_cen[i_halo], cosmo)
        c_cen = mz2c(m200c_cen[i_halo], redshift_cen[i_halo], cosmo)
        for j in 1:n_sh_bar_result[i_halo]
            # i is central halo index, j is index of satellite within each halo
            i_sh = cumush[i_halo]+j # index of satellite in satellite arrays
            parent_sh[i_sh] = i_halo

            log_msat_inner = max(log(rand(T)), T(-7.0))
            r_sh = r_cen * interp.c_lnm2r(
                c_cen, log_msat_inner) * T(200.0^(-1.0/3.0))
            m200c_sh[i_sh] = (interp.muofn(rand(T) * n_sh_bar[i_halo])
                * m200c_cen[i_halo])

            phi = random_phi(T)
            theta = random_theta(T)
            x_sh = halo_pos[1,i_halo] + r_sh * sin(theta) * cos(phi)
            y_sh = halo_pos[2,i_halo] + r_sh * sin(theta) * sin(phi)
            z_sh = halo_pos[3,i_halo] + r_sh * cos(theta)
            dist_sh[i_sh] = sqrt(x_sh^2 + y_sh^2 + z_sh^2)
            redshift_sh[i_sh] = interp.r2z(dist_sh[i_sh])
            theta_sh[i_sh], phi_sh[i_sh] = Healpix.vec2ang(x_sh, y_sh, z_sh)
            lrg_sh[i_sh] = false
            hp_ind_sh[i_sh] = Healpix.vec2pixRing(
                Healpix_res, x_sh, y_sh, z_sh)
        end
        # now how many satellite LRGs do we actually have?
        n_sat_bar = interp.hod_zhengsat(log(m200c_cen[i_halo]))
        n_sat_bar_result = rand(Distributions.Poisson(Float64.(n_sat_bar)))
        sat_perm_sh = sortperm(m200c_sh[cumush[i_halo]+1:cumush[i_halo]+n_sh_bar_result[i_halo]],rev=true)
        for j in 1:n_sat_bar_result
            if j > n_sh_bar_result[i_halo]
                break;
            else
                lrg_sh[cumush[i_halo]+sat_perm_sh[j]] = lrg_cen[i_halo]
            end
        end
    end
end


"""
    generate_sources(model, cosmo, halo_pos_inp, halo_mass_inp; verbose=true)

Produce a source catalog from a model and halo catalog. This converts the
halo arrays into the type specified by `model`.

# Arguments:
- `model::AbstractCIBModel{T}`: source model parameters
- `cosmo::Cosmology.FlatLCDM{T}`: background cosmology
- `Healpix_res::Resolution`: Healpix map resolution
- `halo_pos_inp::AbstractArray{TH,2}`: halo positions with dims (3, nhalos)
- `halo_mass_inp::AbstractArray{TH,1}`: halo masses

# Keywords
- `verbose::Bool=true`: print out progress details
"""
function generate_sources(
        model::AbstractLRGModel{T}, cosmo::Cosmology.FlatLCDM{T},
        halo_pos_inp::AbstractArray{TH,2}, halo_mass_inp::AbstractArray{TH,1},
        halo_quenched::AbstractArray{Bool,1};
        verbose=true) where {T, TH}

    # make sure halo inputs are the CIB type
    halo_pos = convert(Array{T,2}, halo_pos_inp)
    halo_mass = convert(Array{T,1}, halo_mass_inp)

    # set up basics
    N_halos = size(halo_mass, 1)
    interp = get_interpolators( model, cosmo, minimum(halo_mass)*T(0.5), maximum(halo_mass))
    res = Resolution(model.nside)

    verbose && println("Allocating for $(N_halos) centrals.")
    hp_ind_cen = Array{Int64}(undef, N_halos)  # healpix index of halo
    lrg_cen = BitArray(undef, N_halos)  # central has or doesn't have LRG
    redshift_cen = Array{T}(undef, N_halos)
    m200c_cen = Array{T}(undef, N_halos)
    theta_cen = Array{T}(undef, N_halos)
    phi_cen = Array{T}(undef, N_halos)
    dist_cen = Array{T}(undef, N_halos)
    n_sh_bar = Array{T}(undef, N_halos)
    n_sh_bar_result = Array{Int32}(undef, N_halos)

    # STEP 1: compute central properties -----------------------------------
    verbose && println("Processing centrals on $(Threads.nthreads()) threads.")
    process_centrals!(model, cosmo, res,
        interp=interp, hp_ind_cen=hp_ind_cen, dist_cen=dist_cen,
        redshift_cen=redshift_cen, theta_cen=theta_cen, phi_cen=phi_cen, 
        m200c_cen=m200c_cen, lrg_cen=lrg_cen, n_sh_bar=n_sh_bar,
        n_sh_bar_result=n_sh_bar_result,
        halo_pos=halo_pos, halo_mass=halo_mass, halo_quenched=halo_quenched)

    # STEP 2: Generate subhalo arrays -----------------------------
    cumush = generate_subhalo_offsets(n_sh_bar_result)
    total_n_sh = cumush[end]
    hp_ind_sh = Array{Int64}(undef, total_n_sh)  # healpix index of halo
    lrg_sh = BitArray(undef, total_n_sh)  # subhalo has or doesn't have LRG
    redshift_sh = Array{T}(undef, total_n_sh)
    parent_sh = Array{Int64}(undef, total_n_sh)
    m200c_sh = Array{T}(undef, total_n_sh)
    theta_sh = Array{T}(undef, total_n_sh)
    phi_sh = Array{T}(undef, total_n_sh)
    dist_sh = Array{T}(undef, total_n_sh)

    # STEP 3: compute subhalo properties -----------------------------------
    verbose && println("Processing $(total_n_sh) subhalos.")
    process_sats!(model, cosmo, res,
        interp=interp, hp_ind_sh=hp_ind_sh, dist_sh=dist_sh, parent_sh=parent_sh,
        redshift_sh=redshift_sh, theta_sh=theta_sh, phi_sh=phi_sh,
        m200c_sh=m200c_sh, lrg_sh=lrg_sh, cumush=cumush, lrg_cen=lrg_cen,
        m200c_cen=m200c_cen, halo_pos=halo_pos, redshift_cen=redshift_cen,
        n_sh_bar=n_sh_bar, n_sh_bar_result=n_sh_bar_result)
    
    return (
        hp_ind_cen=hp_ind_cen, lrg_cen=lrg_cen, m200c_cen=m200c_cen,
        redshift_cen=redshift_cen, theta_cen=theta_cen, phi_cen=phi_cen, dist_cen=dist_cen,
        hp_ind_sat=hp_ind_sh, lrg_sat=lrg_sh, m200c_sat=m200c_sh, parent_sat=parent_sh,
        redshift_sat=redshift_sh, theta_sat=theta_sh, phi_sat=phi_sh, dist_sat=dist_sh,
        N_cen=N_halos, N_sat=total_n_sh
    )
end

function generate_sources(
        model::AbstractLRGModel{T}, cosmo::Cosmology.FlatLCDM{T},
        halo_pos_inp::AbstractArray{TH,2}, halo_mass_inp::AbstractArray{TH,1};
        verbose=true) where {T, TH}
    dummy_quench = BitArray(undef, length(halo_mass_inp))
    fill!(dummy_quench,true)
    return generate_sources(model,cosmo,halo_pos_inp,halo_mass_inp,dummy_quench;
                             verbose=verbose)
end

"""
    paint!(result_map, model, sources, min_redshift, max_redshift)

'Paint' the LRG catalog onto a map.

# Arguments:
- `result_map::HealpixMap{T_map, RingOrder}`: Healpix map to paint
- `model::AbstractCIBModel{T}`: source model parameters
- `sources`: NamedTuple containing source information from generate_sources
- `fluxes_cen::AbstractArray`: buffer for writing fluxes of centrals
- `fluxes_sat::AbstractArray`: buffer for writing fluxes of satellites
"""
function paint!(result_map::HealpixMap{T_map, RingOrder},
        model::AbstractLRGModel{T}, sources,
        min_redshift::T, max_redshift::T) where {T_map, T}

    pixel_array = result_map.pixels
    fill!(pixel_array, zero(T))  # prepare the frequency map

    # process centrals for this frequency
    Threads.@threads for i in 1:sources.N_cen
        if (sources.lrg_cen[i]) && (sources.redshift_cen[i] < max_redshift) && (sources.redshift_cen[i] > min_redshift)
            pixel_array[sources.hp_ind_cen[i]] += 1
        end
    end

    # process satellites for this frequency
    Threads.@threads for i in 1:sources.N_sat
        if (sources.lrg_sat[i]) && (sources.redshift_sat[i] < max_redshift) && (sources.redshift_sat[i] > min_redshift)
            pixel_array[sources.hp_ind_sat[i]] += 1
        end
    end
end
