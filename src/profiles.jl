

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
compute_θmin(model::AbstractInterpolatorProfile) = exp(first(first(model.ranges)))
compute_θmin(::AbstractProfile{T}) where T = eps(T) 


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


function profile_paint_generic!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}},
                        workspace::CarClenshawCurtisProfileWorkspace, model, α₀, δ₀, 
                        z, Mh, θmax, normalization=1) where T

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
                             normalization * model(θ, z, Mh),
                             zero(T))
        end
    end
end

# fall back to generic profile painter if no specialized painter is defined for the model
function profile_paint!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, 
                        workspace::CarClenshawCurtisProfileWorkspace, model, 
                        α₀, δ₀, z, Mh, θmax, normalization=1) where T
    profile_paint_generic!(m, model, workspace, α₀, δ₀, z, Mh, θmax, normalization)
end


function profile_paint_generic!(m::Enmap{T, 2, Matrix{T}, Gnomonic{T}}, 
                                workspace::GnomonicProfileWorkspace, model, 
                                α₀, δ₀, z, Mh, θmax, normalization=1) where T

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
                             normalization * model(θ, z, Mh),
                             zero(T))
        end
    end
end

# fall back to generic profile painter if no specialized painter is defined for the model
function profile_paint!(m::Enmap{T, 2, Matrix{T}, Gnomonic{T}}, 
                        workspace::GnomonicProfileWorkspace, model, α₀, δ₀, 
                        z, Mh, θmax, normalization=1) where T
    profile_paint_generic!(m, model, workspace, α₀, δ₀, z, Mh, θmax, normalization)
end


function profile_paint_generic!(m::HealpixMap{T, RingOrder}, w::HealpixProfileWorkspace, 
        model, α₀, δ₀, z, Mh, θmax, normalization=1) where T
    ϕ₀ = α₀
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    θmin = max(compute_θmin(model), w.θmin)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(clamp(1 - d² / 2, -one(T), one(T)))
        θ = max(θmin, θ)  # clamp to minimum θ
        m.pixels[ir] += ifelse(θ < θmax, 
                                    normalization * model(θ, z, Mh),
                                    zero(T))
    end
end

# fall back to generic profile painter if no specialized painter is defined for the model
function profile_paint!(m::HealpixMap{T, RingOrder}, w::HealpixProfileWorkspace, model, 
                        α₀, δ₀, z, Mh, θmax, normalization=1) where T
    profile_paint_generic!(m, w, model, α₀, δ₀, z, Mh, θmax, normalization)
end


# paint the the sources in the given range
function paintrange!(m, model, workspace, αs, δs, masses, redshifts, irange::AbstractUnitRange)
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        Mh = masses[i]
        z = redshifts[i]
        θmax_ = compute_θmax(model, Mh * XGPaint.M_sun, z)
        profile_paint!(m, workspace, model, α₀, δ₀, z, Mh, θmax_)
    end
end

# for healpix pixelizations, a buffer is currently required for each thread
function paint!(m::HealpixMap{T, RingOrder}, ws::Vector{W}, model::AbstractProfile, 
        αs, δs, masses, redshifts) where {T, W <: HealpixProfileWorkspace}
    
    fill(m, zero(T))
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize)

    if N_sources < 2Threads.nthreads()  # don't thread if there are not many sources
        return paintrange!(m, first(ws), model, αs, δs, redshifts, masses, 1:N_sources)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paintrange!(m, ws[i], model, αs, δs, redshifts, masses, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paintrange!(m, ws[i], model, αs, δs, redshifts, masses, i1:i2)
    end
end

# staggered threading for safety
function paint!(m, workspace, model, αs, δs, redshifts, masses)
    fill!(m, 0)
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);

    if N_sources < 2Threads.nthreads()  # don't thread if there are not many sources
        return paintrange!(m, workspace, model, αs, δs, redshifts, masses, 1:N_sources)
    end
    
    Threads.@threads :static for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paintrange!(m, workspace, model, αs, δs, redshifts, masses, i1:i2)
    end

    Threads.@threads :static for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paintrange!(m, workspace, model, αs, δs, redshifts, masses, i1:i2)
    end
end
