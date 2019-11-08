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
    a_0::T   = 0.0
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
    II_L_min::T = 2.5e23

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
for each thread. If none are specified, this function will call
`XGPaint.get_thread_RNG()` to create such an array, which takes about 1 second.
"""
function hod_sehgal!(n_I_result, n_II_result,
        halo_mass, redshift, model::Radio_Sehgal2009{T};
        rng_array::Array{<:Random.AbstractRNG,1}=nothing) where T
    N_halos = size(halo_mass,1)

    if rng_array == nothing
        rng_array = get_thread_RNG()
    end
    # compute poisson mean, then draw it
    Threads.@threads for i = 1:N_halos
        I_HON = model.I_N_0 * (halo_mass[i] / model.I_M_0)^model.I_α
        I_HON *= XGPaint.FR_I_redshift_evolution(redshift[i], model)
        n_I_result[i] = rand(rng_array[Threads.threadid()],
            Distributions.Poisson(Float64.(I_HON)))
        II_HON = model.II_N_0 * (halo_mass[i] / model.II_M_0)^model.II_α
        II_HON *= XGPaint.FR_II_redshift_evolution(redshift[i], model)
        n_II_result[i] = rand(rng_array[Threads.threadid()],
            Distributions.Poisson(Float64.(II_HON)))
    end
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
        return Lb * (one(T) + (m * (-one(T) + (Lmin/Lb)^n + n * x))/n)^(one(T)/m)
    end
end

"""
Fills the result array with draws from the luminosity function.
"""
function sehgal_LF!(result::Array{T,1}, m, n, Lb::T, Lmin,
    rng_array::Array{<:Random.AbstractRNG,1}) where T

    if rng_array == nothing
        rng_array = get_thread_RNG()
    end

    if n ≈ 0.0
        Threads.@threads for i = 1:size(result,1)
            result[i] = sehgal_LFn0(
                rand(rng_array[Threads.threadid()]),
                Float64(m), Float64(Lb), Float64(Lmin))
        end
    else
        for i = 1:size(result,1)
            result[i] = sehgal_LF(
                rand(rng_array[Threads.threadid()]),
                Float64(m), Float64(n), Float64(Lb), Float64(Lmin))
        end
    end
end


B(cosθ, β) = ( (1.0-β*cosθ)^(-2.0) + (1.0+β*cosθ)^(-2.0) ) / 2.0

function get_core_lobe_lum(L_beam, ν_Hz, R_int, γ, a_1, a_2, a_3, cosθ)
    β = sqrt(1.0-γ^(-2.0))
    R_obs = R_int * B(cosθ, β)
    L_int = L_beam * (1.0 + R_int) / (1.0 + R_obs)
    L_l_int = L_int / (1.0 + R_int)
    L_c_beam = R_obs * L_l_int

    # now put in the frequency dependence. a_0 = 0.0
    lν = log10(ν_Hz / 1e9 ) # normed to 151 MHz
    f_core =  10 ^ ( a_1 * lν + a_2 * lν^2 + a_3 * lν^3 )

    lν_norm = log10(151e6 / 1e9 ) # normed to 151 MHz
    norm_core =  10 ^ ( a_1 * lν_norm + a_2 * lν_norm^2 + a_3 * lν_norm^3 )

    f_lobe = (ν_Hz / 151e6 ).^(-0.8)
    return L_c_beam * f_core / norm_core, L_l_int * f_lobe
end


function draw_coeff!(a_coeff, model::Radio_Sehgal2009)
    rand!(model.a_1_dist, @view a_coeff[1,:])
    rand!(model.a_2_dist, @view a_coeff[2,:])
    rand!(model.a_3_dist, @view a_coeff[3,:])
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
    rng_array = get_thread_RNG()

    verbose && println("Culling halos below mass $(model.min_mass).")
    mass_cut = halo_mass .> model.min_mass
    halo_pos = halo_pos[:, mass_cut]
    halo_mass = halo_mass[mass_cut]
    N_halos = size(halo_mass, 1)

    verbose && println("Allocating for $(N_halos) halos.")
    hp_ind = Array{Int64}(undef, N_halos)  # healpix index of halo
    redshift = Array{T}(undef, N_halos)
    dist = Array{T}(undef, N_halos)
    r2z = XGPaint.build_r2z_interpolator(
        model.min_redshift, model.max_redshift, cosmo)

    # TODO: make this a utility function
    Threads.@threads for i in 1:N_halos
        dist[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
        redshift[i] = r2z(dist[i])
        hp_ind[i] = Healpix.vec2pixRing(res,
            halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
    end

    verbose && println("Populating HOD.")
    n_I = Array{Int32}(undef, N_halos)
    n_II = Array{Int32}(undef, N_halos)
    hod_sehgal!(n_I, n_II, halo_mass, redshift, model, rng_array=rng_array)

    # set up indices for sources
    cumsat_I = cumsum(n_I)
    prepend!(cumsat_I, 0)
    cumsat_II = cumsum(n_II)
    prepend!(cumsat_II, 0)

    a_coeff_I = Array{T, 2}(undef, 3, cumsat_I[end])
    a_coeff_II = Array{T, 2}(undef, 3, cumsat_II[end])
    draw_coeff!(a_coeff_I, model)
    draw_coeff!(a_coeff_II, model)

    verbose && println("Drawing from luminosity function.")
    # draw luminosities for both populations
    L_I_151 = Array{T, 1}(undef, cumsat_I[end])
    L_II_151 = Array{T, 1}(undef, cumsat_II[end])
    sehgal_LF!(L_I_151, model.I_m, model.I_n, model.I_L_b,
        model.I_L_min, rng_array)
    sehgal_LF!(L_II_151, model.II_m, model.II_n, model.II_L_b,
        model.II_L_min, rng_array)

    return (redshift=redshift, halo_mass=halo_mass, halo_pos=halo_pos,
        dist=dist, n_I=n_I, n_II=n_II,
        a_coeff_I=a_coeff_I, a_coeff_II=a_coeff_II,
        L_I_151=L_I_151, L_II_151=L_II_151
        )
end


export Radio_Sehgal2009
