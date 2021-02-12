"""Functions for computing radio models."""

abstract type AbstractRadioModel{T<:Real} <: AbstractForegroundModel end

"""
    Radio_Sehgal2009{T}(model parameters...)

Define CIB model parameters. Defaults are from Viero et al. 2013.

```@example
model = CIBModel{Float32}(a_0=0.4)
```
"""
Base.@kwdef struct Radio_Sehgal2009{T<:Real} <: AbstractRadioModel{T}
    nside::Int64    = 4096
    min_redshift::T = 0.0
    max_redshift::T = 4.5
    min_mass::T     = 1e13
    box_size::T     = 40000

    # these coefficients are shared for type I and II
    a_0::T   = -1.0
    a_1_dist::Distributions.Uniform{T} = Distributions.Uniform(T(-0.12), T(0.07))
    a_2_dist::Distributions.Uniform{T} = Distributions.Uniform(T(-0.34), T(0.99))
    a_3_dist::Distributions.Uniform{T} = Distributions.Uniform(T(-0.75), T(-0.23))

    I_R_int::T = 10^(-2.6)
    I_γ::T = 6.0
    I_N_0::T = 1.0
    I_M_0::T = 4e13
    I_α::T = 0.1
    I_L_b::T = 10^(24.0)
    I_m::T = -1.55
    I_n::T = 0.0
    I_δ::T = 3.0
    I_z_p::T = 0.8

    II_R_int::T = 10^(-2.8)
    II_γ::T = 8.0
    II_N_0::T = 0.015
    II_M_0::T = 3e15
    II_α::T = 0.1
    II_L_b::T = 10^(27.5)
    II_m::T = -1.6
    II_n::T = -0.65
    II_z_p::T = 1.3
    II_σ_l::T = 0.4
    II_σ_h::T = 0.73

    I_L_min::T = 2.5e23
    II_L_min::T = 4e23

    # physical constants
    phys_h::T     = 6.62606957e-27      # erg.s
    phys_c::T     = 3e+10               # cm/s
    phys_k::T     = 1.3806488e-16       # erg/K
    phys_Msun::T  = 2e33                # g
    phys_Mpc::T   = 3.086e24            # cm
end

pdf_norm1(x, μ, σ) = exp( -(x-μ)^2 / (2 * σ^2) );

function FR_I_redshift_evolution(z, model)
    if z < model.I_z_p
        return (1+z)^model.I_δ
    else
        return (1+model.I_z_p)^model.I_δ
    end
end

function FR_II_redshift_evolution(z, model)
    norm = pdf_norm1(0, model.II_z_p, model.II_σ_l)
    if z < model.II_z_p
        return pdf_norm1(z, model.II_z_p, model.II_σ_l) / norm
    else
        return pdf_norm1(z, model.II_z_p, model.II_σ_h) / norm
    end
end

"""
Populate halos with radio sources according to the HOD in Sehgal et al. 2009.

The optional rng parameter provides an array of random number generators, one
for each thread.
"""
function hod_sehgal(
        halo_mass, redshift, model::Radio_Sehgal2009{T}) where T

    N_halos = size(halo_mass,1)
    nsources_I = Array{Int32}(undef, N_halos)
    nsources_II = Array{Int32}(undef, N_halos)

    # compute poisson mean, then draw it
    Threads.@threads for i = 1:N_halos
        I_HON = model.I_N_0 * (halo_mass[i] / model.I_M_0)^model.I_α
        I_HON *= XGPaint.FR_I_redshift_evolution(redshift[i], model)
        nsources_I[i] = rand(Distributions.Poisson(Float64.(I_HON)))
        II_HON = model.II_N_0 * (halo_mass[i] / model.II_M_0)^model.II_α
        II_HON *= XGPaint.FR_II_redshift_evolution(redshift[i], model)
        nsources_II[i] = rand(Distributions.Poisson(Float64.(II_HON)))
    end

    return nsources_I, nsources_II
end

function sehgal_LFn0(rand_unit::T, m, Lb, Lmin) where T <: Real
    crossover = log( Lb/Lmin )
    x = rand_unit * (-(1/m)+log(Lb/Lmin))
    if x < crossover
        return exp(x) * Lmin
    else
        return Lb * (1+m * x + m * log(Lmin/Lb))^(1/m)
    end
end

function sehgal_LF(rand_unit::T, m, n, Lb, Lmin) where T
    crossover = (T(1) - Lb^(-n) * Lmin^n)/n
    x = (rand_unit * (m-Lb^(-n) * Lmin^n * m-n)/(m * n))
    if x < crossover
        return (Lmin^n + Lb^n * n * x)^(1/n)
    else
        return Lb * (1 + (m * (-1 + (Lmin/Lb)^n + n * x))/n)^(1/m)
    end
end

"""
Fills the result array with draws from the luminosity function.
"""
function sehgal_LF!(result::Array{T,1}, m, n, Lb::T, Lmin,
    rand_array::Array{T,1}) where T

    if n ≈ 0.0
        Threads.@threads for i = 1:size(result,1)
            result[i] = sehgal_LFn0(
                rand_array[i],
                Float64(m), Float64(Lb), Float64(Lmin))
        end
    else
        Threads.@threads for i = 1:size(result,1)
            result[i] = sehgal_LF(
                rand_array[i],
                Float64(m), Float64(n), Float64(Lb), Float64(Lmin))
        end
    end
end


B(cosθ, β) = ( (1-β*cosθ)^(-2) + (1+β*cosθ)^(-2) ) / 2

function get_core_lobe_lum(L_beam::T, ν_Hz, R_int, γ, a_0, a_1, a_2, a_3, cosθ) where T
    β = sqrt(1-γ^(-2))
    R_obs = R_int * B(cosθ, β)
    L_int = L_beam * (1 + R_int) / (1 + R_obs)
    L_l_int = L_int / (1 + R_int)
    L_c_beam = R_obs * L_l_int

    # now put in the frequency dependence. a_0 = 0.0
    lν = log10(ν_Hz / T(1e9) ) # normed to 151 MHz
    f_core =  T(10) ^ ( a_1 * lν + a_2 * lν^2 + a_3 * lν^3 ) # a_0 +

    lν_norm = T(log10(151e6 / 1e9 )) # normed to 151 MHz
    norm_core =  10 ^ ( a_1 * lν_norm + a_2 * lν_norm^2 + a_3 * lν_norm^3 )

    f_lobe = (ν_Hz / T(151e6) ).^T(-0.8)
    return T(L_c_beam * f_core * 0.3f0 / norm_core), T(L_l_int * f_lobe)
end


function draw_coeff!(a_coeff, model::Radio_Sehgal2009)
    a1_task = Threads.@spawn rand!(model.a_1_dist, @view a_coeff[1,:])
    a2_task = Threads.@spawn rand!(model.a_2_dist, @view a_coeff[2,:])
    a3_task = Threads.@spawn rand!(model.a_3_dist, @view a_coeff[3,:])

    wait(a1_task)
    wait(a2_task)
    wait(a3_task)
end

"""
Produce a source catalog from a model and halo catalog.
"""
function generate_sources(
        # model parameters
        model::AbstractRadioModel, cosmo::Cosmology.FlatLCDM{T},
        # halo arrays
        halo_pos::Array{T,2}, halo_mass::Array{T,1};
        verbose=true) where T

    res = Resolution(model.nside)

    verbose && println("Culling halos below mass $(model.min_mass).")
    mass_cut = halo_mass .> model.min_mass
    halo_pos = halo_pos[:, mass_cut]
    halo_mass = halo_mass[mass_cut]
    N_halos = size(halo_mass, 1)

    verbose && println("Allocating for $(N_halos) halos.")
    dist, redshift, hp_ind = get_basic_halo_properties(halo_pos, model, cosmo, res)
    θ, ϕ = get_angles(halo_pos)

    verbose && println("Populating HOD.")
    nsources_I, nsources_II = hod_sehgal(
        halo_mass, redshift, model)

     # total number of sources
    total_n_I, total_n_II = sum(nsources_I), sum(nsources_II)

    # the last element of cumsat_I is the total number of type I sources
    verbose && println("Drawing spectral coefficients.")
    a_coeff_I = Array{T, 2}(undef, 3, total_n_I)
    a_coeff_II = Array{T, 2}(undef, 3, total_n_II)
    draw_coeff!(a_coeff_I, model)
    draw_coeff!(a_coeff_II, model)

    verbose && println("Drawing from luminosity function.")
    # generate random numbers for LF
    rand_buffer_I = Array{T, 1}(undef, total_n_I)
    rand_buffer_II = Array{T, 1}(undef, total_n_II)
    threaded_rand!(rand_buffer_I)
    threaded_rand!(rand_buffer_II)

    # draw luminosities for both populations
    L_I_151 = Array{T, 1}(undef, total_n_I)
    L_II_151 = Array{T, 1}(undef, total_n_II)
    sehgal_LF!(L_I_151, model.I_m, model.I_n, model.I_L_b,
        model.I_L_min, rand_buffer_I)
    sehgal_LF!(L_II_151, model.II_m, model.II_n, model.II_L_b,
        model.II_L_min, rand_buffer_II)

    verbose && println("Drawing for impact parameter.")
    # reuse rand buffers for the impact parameters
    threaded_rand!(rand_buffer_I)
    threaded_rand!(rand_buffer_II)
    cosθ_I = rand_buffer_I
    cosθ_II = rand_buffer_II

    return (
        redshift=redshift, halo_mass=halo_mass, halo_pos=halo_pos,
        halo_hp_ind=hp_ind, θ=θ, ϕ=ϕ,
        dist=dist, nsources_I=nsources_I, nsources_II=nsources_II,
        a_coeff_I=a_coeff_I, a_coeff_II=a_coeff_II,
        L_I_151=L_I_151, L_II_151=L_II_151,
        cosθ_I=cosθ_I, cosθ_II=cosθ_II,
        total_n_I=total_n_I, total_n_II=total_n_II
    )
end


function paint!(result_map::Map{T_map,RingOrder},
                nu_obs::T, model::AbstractRadioModel, sources;
                return_fluxes=false) where {T_map, T}

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
            L_c, L_l = XGPaint.get_core_lobe_lum(
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
            L_c, L_l = XGPaint.get_core_lobe_lum(
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


export Radio_Sehgal2009, paint!
