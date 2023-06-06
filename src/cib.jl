"""Functions for computing CIB models."""

abstract type AbstractCIBModel{T<:Real} <: AbstractForegroundModel end

"""
    CIB_Planck2013{T}(; kwargs...)

Define CIB model parameters. Defaults are from Viero et al. 2013. All numbers
not typed are converted to type T. This model has the following parameters and default values:

* `nside::Int64 = 4096`
* `hod::String = "shang"`
* `Inu_norm = 0.3180384`
* `min_redshift = 0.0`
* `max_redshift = 5.0`
* `min_mass = 1e12`
* `box_size = 40000`
* `shang_zplat = 2.0`
* `shang_Td = 20.7`
* `shang_betan = 1.6`
* `shang_eta = 2.4`
* `shang_alpha = 0.2`
* `shang_Mpeak = 10^12.3`
* `shang_sigmaM = 0.3`
* `shang_Msmin = 1e11`
* `shang_Mmin = 1e10`
* `shang_I0 = 46`
* `jiang_gamma_1 = 0.13`
* `jiang_alpha_1 = -0.83`
* `jiang_gamma_2 = 1.33`
* `jiang_alpha_2 = -0.02`
* `jiang_beta_2 = 5.67`
* `jiang_zeta = 1.19`

"""
@with_kw struct CIB_Planck2013{T<:Real} <: AbstractCIBModel{T} @deftype T
    nside::Int64    = 4096
    hod::String     = "shang"
    Inu_norm     = 0.3180384
    min_redshift = 0.0
    max_redshift = 5.0
    min_mass     = 1e12
    box_size     = 40000

    # shang HOD
    shang_zplat  = 2.0
    shang_Td     = 20.7
    shang_beta   = 1.6
    shang_eta    = 2.4
    shang_alpha  = 0.2
    shang_Mpeak  = 10^12.3
    shang_sigmaM = 0.3
    shang_Msmin  = 1e11
    shang_Mmin   = 1e10
    shang_I0     = 46

    # jiang
    jiang_gamma_1    = 0.13
    jiang_alpha_1    = -0.83
    jiang_gamma_2    = 1.33
    jiang_alpha_2    = -0.02
    jiang_beta_2     = 5.67
    jiang_zeta       = 1.19
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
        r2z = build_r2z_interpolator(
            model.min_redshift, model.max_redshift, cosmo),
        hod_shang = build_shang_interpolator(
            log(min_halo_mass), log(max_halo_mass), model),
        c_lnm2r = build_c_lnm2r_interpolator(),
        sigma_sat = build_sigma_sat_ln_interpolator(
            log(max_halo_mass), model),
        muofn = build_muofn_interpolator(model)
    )
end



# Fill up arrays with information related to CIB central sources.
function process_centrals!(
    model::AbstractCIBModel{T}, cosmo::Cosmology.FlatLCDM{T}, Healpix_res::Resolution;
    interp, hp_ind_cen, dist_cen, redshift_cen, theta_cen, phi_cen,
    lum_cen, n_sat_bar, n_sat_bar_result,
    halo_pos, halo_mass) where T

    N_halos = size(halo_mass, 1)

    Threads.@threads for i = 1:N_halos
        # location information for centrals
        hp_ind_cen[i] = Healpix.vec2pixRing(Healpix_res,
            halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
        theta_cen[i], phi_cen[i] = Healpix.vec2ang(halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
        dist_cen[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
        redshift_cen[i] = interp.r2z(dist_cen[i])

        # compute HOD
        n_sat_bar[i] = interp.hod_shang(log(halo_mass[i]))
        n_sat_bar_result[i] = rand(Distributions.Poisson(Float64.(n_sat_bar[i])))

        # get central luminosity
        lum_cen[i] = sigma_cen(halo_mass[i], model) * shang_z_evo(
            redshift_cen[i], model)
    end
end


# Fill up arrays with information related to CIB satellites.
function process_sats!(
        model::AbstractCIBModel{T}, cosmo::Cosmology.FlatLCDM{T},
        Healpix_res::Resolution;
        interp, hp_ind_sat, dist_sat, redshift_sat, theta_sat, phi_sat,
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
            theta_sat[i_sat], phi_sat[i_sat] = Healpix.vec2ang(x_sat, y_sat, z_sat)

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
        model::AbstractCIBModel{T}, cosmo::Cosmology.FlatLCDM{T},
        halo_pos_inp::AbstractArray{TH,2}, halo_mass_inp::AbstractArray{TH,1};
        verbose=true) where {T, TH}

    # make sure halo inputs are the CIB type
    halo_pos = convert(Array{T,2}, halo_pos_inp)
    halo_mass = convert(Array{T,1}, halo_mass_inp)

    # set up basics
    N_halos = size(halo_mass, 1)
    interp = get_interpolators( model, cosmo, minimum(halo_mass), maximum(halo_mass))
    res = Resolution(model.nside)

    verbose && println("Allocating for $(N_halos) centrals.")
    hp_ind_cen = Array{Int64}(undef, N_halos)  # healpix index of halo
    lum_cen = Array{T}(undef, N_halos)  # Lum of central w/o ν-dependence
    redshift_cen = Array{T}(undef, N_halos)
    theta_cen = Array{T}(undef, N_halos)
    phi_cen = Array{T}(undef, N_halos)
    dist_cen = Array{T}(undef, N_halos)
    n_sat_bar = Array{T}(undef, N_halos)
    n_sat_bar_result = Array{Int32}(undef, N_halos)

    # STEP 1: compute central properties -----------------------------------
    verbose && println("Processing centrals on $(Threads.nthreads()) threads.")
    process_centrals!(model, cosmo, res,
        interp=interp, hp_ind_cen=hp_ind_cen, dist_cen=dist_cen,
        redshift_cen=redshift_cen, theta_cen=theta_cen, phi_cen=phi_cen, 
        lum_cen=lum_cen, n_sat_bar=n_sat_bar,
        n_sat_bar_result=n_sat_bar_result,
        halo_pos=halo_pos, halo_mass=halo_mass)

    # STEP 2: Generate satellite arrays -----------------------------
    cumsat = generate_subhalo_offsets(n_sat_bar_result)
    total_n_sat = cumsat[end]
    hp_ind_sat = Array{Int64}(undef, total_n_sat)  # healpix index of halo
    lum_sat = Array{T}(undef, total_n_sat)  # Lum of central w/o ν-dependence
    redshift_sat = Array{T}(undef, total_n_sat)
    theta_sat = Array{T}(undef, total_n_sat)
    phi_sat = Array{T}(undef, total_n_sat)
    dist_sat = Array{T}(undef, total_n_sat)

    # STEP 3: compute satellite properties -----------------------------------
    verbose && println("Processing $(total_n_sat) satellites.")
    process_sats!(model, cosmo, res,
        interp=interp, hp_ind_sat=hp_ind_sat, dist_sat=dist_sat,
        redshift_sat=redshift_sat, theta_sat=theta_sat, phi_sat=phi_sat,
        lum_sat=lum_sat, cumsat=cumsat,
        halo_mass=halo_mass, halo_pos=halo_pos, redshift_cen=redshift_cen,
        n_sat_bar=n_sat_bar, n_sat_bar_result=n_sat_bar_result)
    
    return (
        hp_ind_cen=hp_ind_cen, lum_cen=lum_cen,
        redshift_cen=redshift_cen, theta_cen=theta_cen, phi_cen=phi_cen, dist_cen=dist_cen,
        hp_ind_sat=hp_ind_sat, lum_sat=lum_sat,
        redshift_sat=redshift_sat, theta_sat=theta_sat, phi_sat=phi_sat, dist_sat=dist_sat,
        N_cen=N_halos, N_sat=total_n_sat
    )
end



function fill_fluxes!(nu_obs, model::AbstractCIBModel{T}, sources,
        fluxes_cen::AbstractArray, fluxes_sat::AbstractArray) where T
        
    # process centrals for this frequency
    Threads.@threads for i in 1:sources.N_cen
        nu = (one(T) + sources.redshift_cen[i]) * nu_obs
        fluxes_cen[i] = l2f(
            sources.lum_cen[i] * nu2theta(
                nu, sources.redshift_cen[i], model),
            sources.dist_cen[i], sources.redshift_cen[i])
    end

    # process satellites for this frequency
    Threads.@threads for i in 1:sources.N_sat
        nu = (one(T) + sources.redshift_sat[i]) * nu_obs
        fluxes_sat[i] = l2f(
            sources.lum_sat[i] * nu2theta(
                nu, sources.redshift_sat[i], model),
            sources.dist_sat[i], sources.redshift_sat[i])
    end
end


"""
    paint!(result_map, nu_obs, model, sources, fluxes_cen, fluxes_sat)

Paint a source catalog onto a map, recording the fluxes in
`fluxes_cen` and `fluxes_sat`.

# Arguments:
- `result_map::HealpixMap{T_map, RingOrder}`: Healpix map to paint
- `nu_obs`: frequency in Hz
- `model::AbstractCIBModel{T}`: source model parameters
- `sources`: NamedTuple containing source information from generate_sources
- `fluxes_cen::AbstractArray`: buffer for writing fluxes of centrals
- `fluxes_sat::AbstractArray`: buffer for writing fluxes of satellites
"""
function paint!(result_map::HealpixMap{T_map, RingOrder},
        nu_obs, model::AbstractCIBModel{T}, sources,
        fluxes_cen::AbstractArray, fluxes_sat::AbstractArray) where {T_map, T}

    pixel_array = result_map.pixels
    fill!(pixel_array, zero(T))  # prepare the frequency map

    fill_fluxes!(nu_obs, model, sources, fluxes_cen, fluxes_sat)

    # process centrals for this frequency
    Threads.@threads for i in 1:sources.N_cen
        pixel_array[sources.hp_ind_cen[i]] += fluxes_cen[i]
    end

    # process satellites for this frequency
    Threads.@threads for i in 1:sources.N_sat
        pixel_array[sources.hp_ind_sat[i]] += fluxes_sat[i]
    end

    # divide by healpix pixel size
    per_pixel_steradian = 1 / nside2pixarea(result_map.resolution.nside)
    pixel_array .*= per_pixel_steradian
    return result_map
end

# CAR version
function paint!(result_map::Enmap{TM},
        nu_obs, model::AbstractCIBModel{T}, sources,
        fluxes_cen::AbstractArray, fluxes_sat::AbstractArray) where {TM, T}

    fill!(result_map, zero(TM))  # zero out the frequency map
    fill_fluxes!(nu_obs, model, sources, fluxes_cen, fluxes_sat)  # generate fluxes

    pixsizes = pixareamap(result_map)
    catalog2map!(result_map, fluxes_cen, sources.theta_cen, sources.phi_cen, pixsizes, erase_first=false)
    catalog2map!(result_map, fluxes_sat, sources.theta_sat, sources.phi_sat, pixsizes, erase_first=false)
    return result_map
end





"""
Paint a source catalog onto a map.

This function creates the arrays for you.
"""
function paint!(result_map::HealpixMap{T,RingOrder},
        nu_obs::T, model::AbstractCIBModel, sources) where T

    fluxes_cen = Array{T, 1}(undef, sources.N_cen)
    fluxes_sat = Array{T, 1}(undef, sources.N_sat)

    paint!(result_map, nu_obs, model, sources,
        fluxes_cen, fluxes_sat)

    return fluxes_cen, fluxes_sat
end
