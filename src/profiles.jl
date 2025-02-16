

# RECTANGULAR WORKSPACES

abstract type AbstractProfileWorkspace end

struct CarClenshawCurtisProfileWorkspace{A} <: AbstractProfileWorkspace
    sin_α::A
    cos_α::A
    sin_δ::A
    cos_δ::A
end

function profileworkspace(shape, wcs::CarClenshawCurtis)
    α_map, δ_map = posmap(shape, wcs)
    return CarClenshawCurtisProfileWorkspace(
        sin.(α_map), cos.(α_map), sin.(δ_map), cos.(δ_map))
end

struct GnomonicProfileWorkspace{A} <: AbstractProfileWorkspace
    sin_α::A
    cos_α::A
    sin_δ::A
    cos_δ::A
end

function profileworkspace(shape, wcs::Gnomonic)
    α_map, δ_map = posmap(shape, wcs)
    return GnomonicProfileWorkspace(
        sin.(α_map), cos.(α_map), sin.(δ_map), cos.(δ_map))
end


abstract type AbstractProfile{T} end
abstract type AbstractGNFW{T} <: AbstractProfile{T} end

struct Battaglia16ThermalSZProfile{T,C} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
end

struct BreakModel{T,C} <: AbstractGNFW{T}
    f_b::T
    cosmo::C
    alpha_break::T
    M_break::T
end

function Battaglia16ThermalSZProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    return Battaglia16ThermalSZProfile(f_b, cosmo)
end

function BreakModel(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, alpha_break::T=1.5, M_break::T=2.0*10^14) where {T <: Real}
    #alpha_break = 1.486 from Shivam P paper by Nate's sleuthing
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    return BreakModel(f_b, cosmo, alpha_break, M_break)
end

const ρ_crit_factor = uconvert(u"kg/m^3", 3u"km^2*Mpc^-2*s^-2" / (8π * constants.G))


function ρ_crit(model, z)
    H_z = H(model.cosmo, z)
    return uconvert(u"kg/m^3", 3H_z^2 / (8π * constants.G))
end

function R_Δ(model, M_Δ, z, Δ=200)
    return ∛(M_Δ / (4π/3 * Δ * ρ_crit(model, z)))
end

function angular_size(model::AbstractProfile{T}, physical_size, z) where T
    d_A = angular_diameter_dist(model.cosmo, z)

    # convert both to the same units and strip units for atan
    phys_siz_unitless = T(ustrip(uconvert(unit(d_A), physical_size)))
    d_A_unitless = T(ustrip(d_A))
    return atan(phys_siz_unitless, d_A_unitless)
end

function generalized_nfw(x, xc, α, β, γ)
    x̄ = x / xc
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end

function _generalized_scaled_nfw(x̄, α, β, γ)
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end


function get_params(::AbstractGNFW{T}, M_200, z) where T
	z₁ = z + 1
	m = M_200 / (1e14M_sun)
	P₀ = 18.1 * m^0.154 * z₁^-0.758
	xc = 0.497 * m^-0.00865 * z₁^0.731
	β = 4.35 * m^0.0393 * z₁^0.415
	α = 1
    γ = -0.3
    β = γ - α * β  # Sigurd's conversion from Battaglia to standard NFW
    return (xc=T(xc), α=T(α), β=T(β), γ=T(γ), P₀=T(P₀))
end

# _tsz_y₁(x, _a) = (x*(_a+1))^(1/(_a+1))
# _tsz_x₁(y, _a) = y^(_a+1)/(_a+1)
function _nfw_profile_los_quadrature(x, xc, α, β, γ; zmax=1e5, rtol=eps(), order=9)
    x² = x^2
    scale = 1e9
    integral, err = quadgk(y -> scale * generalized_nfw(√(y^2 + x²), xc, α, β, γ),
                      0.0, zmax, rtol=rtol, order=order)
    return 2integral / scale
end

function dimensionless_P_profile_los(model::Battaglia16ThermalSZProfile{T}, M_200, z, r) where T
    par = get_params(model, M_200, z)
    R_200 = R_Δ(model, M_200, z, 200)
    x = r / angular_size(model, R_200, z)
    return par.P₀ * _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
end

function dimensionless_P_profile_los(model::BreakModel{T}, M_200, z, r) where T
    par = get_params(model, M_200, z)
    R_200 = R_Δ(model, M_200, z, 200)
    x = r / angular_size(model, R_200, z)
    if M_200 < model.M_break * M_sun
        return par.P₀ * (M_200/(model.M_break*M_sun))^model.alpha_break * _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
    else
        return par.P₀ * _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
    end
end

"""Line-of-sight integrated electron pressure"""
P_e_los(model, M_200, z, r) = 0.5176 * P_th_los(model, M_200, z, r)

"""Line-of-sight integrated thermal pressure"""
P_th_los(model, M_200, z, r) = constants.G * M_200 * 200 * ρ_crit(model, z) * 
    model.f_b / 2 * dimensionless_P_profile_los(model, M_200, z, r)

function compton_y(model, M_200, z, r)
    return P_e_los(model, M_200, z, r) * P_e_factor
end

(model::Battaglia16ThermalSZProfile)(r, M_200, z) = compton_y(model, M_200, z, r)

function profile_grid(model::AbstractGNFW{T}; N_z=256, N_logM=256, N_logθ=512, z_min=1e-3, 
        z_max=5.0, logM_min=11, logM_max=15.7, logθ_min=-16.5, logθ_max=2.5) where T

    logθs = LinRange(logθ_min, logθ_max, N_logθ)
    redshifts = LinRange(z_min, z_max, N_z)
    logMs = LinRange(logM_min, logM_max, N_logM)
    return profile_grid(model, logθs, redshifts, logMs)
end

function profile_grid(model::AbstractGNFW{T}, logθs, redshifts, logMs) where T

    N_logθ, N_z, N_logM = length(logθs), length(redshifts), length(logMs)
    A = zeros(T, (N_logθ, N_z, N_logM))

    Threads.@threads :static for im in 1:N_logM
        logM = logMs[im]
        M = 10^(logM) * M_sun
        for (iz, z) in enumerate(redshifts)
            for iθ in 1:N_logθ
                θ = exp(logθs[iθ])
                A[iθ, iz, im] = max(zero(T), model(M, z, θ))
            end
        end
    end

    return logθs, redshifts, logMs, A
end


# get angular size in radians of radius to stop at
function compute_θmax(model::AbstractProfile{T}, M_Δ, z; mult=4) where T
    r = R_Δ(model, M_Δ, z)
    return T(mult * angular_size(model, r, z))
end

# prevent infinities at cusp
compute_θmin(model::AbstractLogInterpolatorProfile) = exp(first(first(model.ranges)))
compute_θmin(::AbstractProfile{T}) where T = eps(T) 

# DEBUGGING ONLY: VERY APPROXIMATE
function websky_m200m_to_m200c(m200m, z, cosmo)
    Ω_m = cosmo.Ω_m
    omz = Ω_m * (1+z)^3 / ( Ω_m * (1+z)^3 + 1 - Ω_m )
    m200c = omz^0.35 * m200m

    return m200c
end

# find maximum radius to integrate to
function build_max_paint_logradius(logθs, redshifts, logMs, 
                              A::AbstractArray{T}; rtol=1e-2) where T
    
    logRs = zeros(T, (size(A)[2:3]))
    N_logM = length(logMs)
    N_logθ = length(logθs)
    dF_r = zeros(N_logθ)
    
    for im in 1:N_logM
        for (iz, z) in enumerate(redshifts)
            s = zero(T)
            for iθ in 1:(N_logθ-1)
                θ₁ = exp(logθs[iθ])
                θ₂ = exp(logθs[iθ+1])
                f₁ = A[iθ, iz, im] * θ₁
                f₂ = A[iθ+1, iz, im] * θ₂
                s += (θ₂ - θ₁) * (f₁ + f₂) / 2

                dF_r[iθ] = s
            end

            threshold = (1-rtol) * s
            for iθ in (N_logθ-1):-1:1
                if dF_r[iθ] < threshold
                    logRs[iz, im] = min(logθs[iθ], log(π))
                    break
                end
            end
            
        end
    end

    return scale(
        Interpolations.interpolate(logRs, BSpline(Cubic(Line(OnGrid())))), 
        redshifts, logMs);
end

abstract type AbstractInterpolatorProfile{T} <: AbstractProfile{T} end


"""
    LogInterpolatorProfile{T, P, I1}

A profile that interpolates over a positive-definite function (θ, z, M_halo), but internally
interpolates over log(θ) and log10(M) using a given interpolator. Evaluation of this profile
is then done by exponentiating the result of the interpolator.

```
    f(θ, z, M) = exp(itp(log(θ), z, log10(M)))
```

This is useful for interpolating over a large range of scales and masses, where the profile
is expected to be smooth in log-log space. It wraps the original model and also the 
interpolator object itself.
"""
struct LogInterpolatorProfile{T, P <: AbstractProfile{T}, I1} <: AbstractInterpolatorProfile{T}
    model::P
    itp::I1
end

# forward the interpolator calls to the wrapped interpolator
@inline (ip::LogInterpolatorProfile)(θ, z, Mh_Msun) = exp(ip.itp(log(θ), z, log10(Mh_Msun)))

Base.show(io::IO, ip::LogInterpolatorProfile{T,P,I1}) where {T,P,I1} = print(
    io, "LogInterpolatorProfile{$(T),\n  $(P),\n  ...} interpolating over size ", size(ip.itp))


"""Helper function to build a (θ, z, Mh) interpolator"""
function build_interpolator(model::AbstractProfile; cache_file::String="", 
                            N_logθ=512, pad=256, overwrite=true, verbose=true)

    if overwrite || (isfile(cache_file) == false)
        verbose && print("Building new interpolator from model.\n")
        rft = RadialFourierTransform(n=N_logθ, pad=pad)
        logθ_min, logθ_max = log(minimum(rft.r)), log(maximum(rft.r))
        prof_logθs, prof_redshift, prof_logMs, prof_y = profile_grid(model; 
            N_logθ=N_logθ, logθ_min=logθ_min, logθ_max=logθ_max)
        if length(cache_file) > 0
            verbose && print("Saving new interpolator to $(cache_file).\n")
            save(cache_file, Dict("prof_logθs"=>prof_logθs, 
                "prof_redshift"=>prof_redshift, "prof_logMs"=>prof_logMs, "prof_y"=>prof_y))
        end
    else
        print("Found cached Battaglia profile model. Loading from disk.\n")
        model_grid = load(cache_file)
        prof_logθs, prof_redshift, prof_logMs, prof_y = model_grid["prof_logθs"], 
            model_grid["prof_redshift"], model_grid["prof_logMs"], model_grid["prof_y"]
    end

    itp = Interpolations.interpolate(log.(prof_y), BSpline(Cubic(Line(OnGrid()))))
    interp_model = scale(itp, prof_logθs, prof_redshift, prof_logMs)
    return LogInterpolatorProfile(model, interp_model)
end


function profile_paint!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, 
                        α₀, δ₀, workspace::CarClenshawCurtisProfileWorkspace, 
                        model, z, Mh, θmax, mult_factor=1) where T

    # get indices of the region to work on
    i1, j1 = sky2pix(m, α₀ - θmax, δ₀ - θmax)
    i2, j2 = sky2pix(m, α₀ + θmax, δ₀ + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))
    θmin = compute_θmin(model)

    x₀ = cos(δ₀) * cos(α₀)
    y₀ = cos(δ₀) * sin(α₀) 
    z₀ = sin(δ₀)

    @inbounds for j in j_start:j_stop
        for i in i_start:i_stop
            x₁ = workspace.cos_δ[i,j] * workspace.cos_α[i,j]
            y₁ = workspace.cos_δ[i,j] * workspace.sin_α[i,j]
            z₁ = workspace.sin_δ[i,j]
            d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
            θ =  acos(clamp(1 - d² / 2, -one(T), one(T)))
            θ = max(θmin, θ)  # clamp to minimum θ
            m[i,j] += ifelse(θ < θmax, 
                             mult_factor * model(θ, z, Mh),
                             zero(T))
        end
    end
end


function profile_paint!(m::Enmap{T, 2, Matrix{T}, Gnomonic{T}}, 
            α₀, δ₀, workspace::GnomonicProfileWorkspace, model, z, Ms, θmax, mult_factor=1) where T

    # get indices of the region to work on
    i1, j1 = sky2pix(m, α₀ - θmax, δ₀ - θmax)
    i2, j2 = sky2pix(m, α₀ + θmax, δ₀ + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))
    θmin = compute_θmin(model)

    x₀ = cos(δ₀) * cos(α₀)
    y₀ = cos(δ₀) * sin(α₀) 
    z₀ = sin(δ₀)

    @inbounds for j in j_start:j_stop
        for i in i_start:i_stop
            x₁ = workspace.cos_δ[i,j] * workspace.cos_α[i,j]
            y₁ = workspace.cos_δ[i,j] * workspace.sin_α[i,j]
            z₁ = workspace.sin_δ[i,j]
            d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
            θ =  acos(clamp(1 - d² / 2, -one(T), one(T)))
            θ = max(θmin, θ)  # clamp to minimum θ
            m[i,j] += ifelse(θ < θmax, 
                             mult_factor * model(θ, z, Mh),
                             zero(T))
        end
    end
end


function profile_paint!(m::HealpixMap{T, RingOrder}, 
            α₀, δ₀, w::HealpixProfileWorkspace, model, z, Mh, θmax, mult_factor=1) where T
    ϕ₀ = α₀
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    θmin = max(compute_θmin(model), w.θmin)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ =  acos(clamp(1 - d² / 2, -one(T), one(T)))
        θ = max(θmin, θ)  # clamp to minimum θ
        m.pixels[ir] += ifelse(θ < θmax, 
                                    mult_factor * model(θ, z, Mh),
                                    zero(T))
    end
end


# for rectangular pixelizations

# multi-halo painting utilities
function paint!(m, model, workspace, 
                masses::AV, redshifts::AV, αs::AV, δs::AV, irange::AbstractUnitRange) where AV
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax_ = compute_θmax(p, mh * XGPaint.M_sun, z)
        profile_paint!(m, α₀, δ₀, workspace, model, z, mh, θmax_)
    end
end


function paint!(m::HealpixMap{T, RingOrder}, p::XGPaint.AbstractProfile, ws::Vector{W}, interp, masses::AV, 
                    redshifts::AV, αs::AV, δs::AV) where {T, W <: HealpixProfileWorkspace, AV}
    m .= 0.0

    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint!(m, p, ws[i], interp, masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint!(m, p, ws[i], interp, masses, redshifts, αs, δs, i1:i2)
    end
end

function paint!(m, p::XGPaint.AbstractProfile, workspace, model, masses::AV, 
                        redshifts::AV, αs::AV, δs::AV)  where AV
    fill!(m, 0)
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);

    if N_sources < 2Threads.nthreads()  # don't thread if there are not many sources
        return paint!(m, model, workspace, masses, redshifts, αs, δs, 1:N_sources)
    end
    
    Threads.@threads :static for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint!(m, model, workspace, masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads :static for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint!(m, model, workspace, masses, redshifts, αs, δs, i1:i2)
    end
end
