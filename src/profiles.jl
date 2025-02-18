

# RECTANGULAR WORKSPACES


# default workspaces are immutable, so just forward the type
wrapserialworkspace(w, tid) = w

struct CarClenshawCurtisProfileWorkspace{T,A<:AbstractArray{T,2}} <: AbstractProfileWorkspace{T}
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


struct GnomonicProfileWorkspace{T,A<:AbstractArray{T,2}} <: AbstractProfileWorkspace{T}
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

Base.show(io::IO, w::AbstractProfileWorkspace) = print(io, "$(typeof(w))")


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
        M = 10^(logM)
        for (iz, z) in enumerate(redshifts)
            for iθ in 1:N_logθ
                θ = exp(logθs[iθ])
                A[iθ, iz, im] = max(zero(T), model(θ, M, z))
            end
        end
    end

    return logθs, redshifts, logMs, A
end



"""
Computes a real-space beam interpolator and a maximum
"""
function realspacegaussbeam(::Type{T}, θ_FWHM::Ti; rtol=1e-24, N_θ::Int=2000) where {T,Ti}
    Nlmax = ceil(Int, log2(8π / θ_FWHM))
    lmax = 2^Nlmax

    b_l = gaussbeam(θ_FWHM, lmax)
    θs = LinRange(zero(θ_FWHM), 5θ_FWHM, N_θ)
    b_θ = XGPaint.bl2beam(b_l, θs)
    atol = b_θ[begin] * rtol
    i_max = findfirst(<(atol), b_θ)

    θs = convert(LinRange{T, Int}, θs[begin:i_max])
    beam_real_interp = cubic_spline_interpolation(
        θs, T.(b_θ[begin:i_max]), extrapolation_bc=zero(T))
    return beam_real_interp, θs
end


# heuristic for sizehint
_approx_discbuffer_size(nside, θmax) = ceil(Int, 1.1 * π * θmax^2 / (nside2pixarea(nside)))

struct HealpixSerialProfileWorkspace{T} <: AbstractProfileWorkspace{T}
    nside::Int
    ringinfo::RingInfo
    disc_buffer::Vector{Int}
    θmax::T
    posmap::Vector{Tuple{T, T, T}}

    HealpixSerialProfileWorkspace{T}(nside::Int, θmax) where T = new{T}(
        nside,
        RingInfo(0, 0, 0, 0.0, true),
        sizehint!(Int[], _approx_discbuffer_size(nside, θmax)),
        T(θmax),
        vectorhealpixmap(T, nside)
    )

    HealpixSerialProfileWorkspace{T}(nside, ringinfo, disc_buffer, θmax, posmap) where T = new{T}(
        nside, ringinfo, disc_buffer, T(θmax), posmap)
end

# default outer constructors for Float64
HealpixSerialProfileWorkspace(nside::Int, θmax) = HealpixSerialProfileWorkspace{Float64}(nside, θmax)

function Base.show(io::IO, ::HealpixSerialProfileWorkspace{T}) where T
    expr = "HealpixSerialProfileWorkspace{$(T)}"
    print(io, expr)
end

struct HealpixProfileWorkspace{T}  <: AbstractProfileWorkspace{T}
    nside::Int
    ringinfo::RingInfo
    disc_buffer::Vector{Vector{Int}}
    θmax::T
    posmap::Vector{Tuple{T, T, T}}

    HealpixProfileWorkspace{T}(nside::Int, θmax, nthreads=Threads.nthreads()) where T = new{T}(
        nside,
        RingInfo(0, 0, 0, 0.0, true),
        Vector{Int}[sizehint!(Int[], _approx_discbuffer_size(nside, θmax)) for i in 1:nthreads],
        T(θmax),
        vectorhealpixmap(T, nside)
    )
end

HealpixProfileWorkspace(nside::Int, θmax) = HealpixProfileWorkspace{Float64}(nside, θmax)

# unlike other workspaces, Healpix workspace needs a mutable buffer per thread
# so asking for a workspace with a thread id will return the appropriate mutable workspace
function wrapserialworkspace(w::HealpixProfileWorkspace{T}, tid) where T
    return HealpixSerialProfileWorkspace{T}(w.nside, w.ringinfo, w.disc_buffer[tid], w.θmax, w.posmap)
end


function realspacebeampaint!(hp_map, w::HealpixSerialProfileWorkspace, realprofile, flux, θ₀, ϕ₀)
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, hp_map.resolution, θ₀, ϕ₀, w.θmax)

    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(1 - d² / 2)
        hp_map.pixels[ir] += flux * realprofile(θ)
    end
end


"""Apply a beam to a profile grid"""
function transform_profile_grid!(y_prof_grid, rft, lbeam)
    rprof = y_prof_grid[:,1,1]
    for i in axes(y_prof_grid, 2)
        for j in axes(y_prof_grid, 3)
            rprof .= y_prof_grid[:,i,j]
            lprof = real2harm(rft, rprof)
            lprof .*= lbeam
            reverse!(lprof)
            rprof′ = harm2real(rft, lprof)
            y_prof_grid[:,i,j] .= rprof′
        end
    end
end

"prune a profile grid for negative values, extrapolate instead"
function cleanup_negatives!(y_prof_grid)
    for i in axes(y_prof_grid, 2)
        for j in axes(y_prof_grid, 3)
            extrapolating = false
            fact = 1.0
            for k in axes(y_prof_grid, 1)
                if y_prof_grid[k,i,j] <= 0
                    extrapolating = true
                    fact = y_prof_grid[k-1,i,j] / y_prof_grid[k-2,i,j]
                end
                if extrapolating
                    y_prof_grid[k,i,j] = max(fact * y_prof_grid[k-1,i,j], nextfloat(0.0))
                end
            end
        end
    end
end



# get angular size in radians of radius to stop at
function compute_θmax(model::AbstractProfile{T}, M_Δ, z; mult=4) where T
    r = R_Δ(model, M_Δ, z)
    return T(mult * angular_size(model, r, z))
end

# prevent infinities at cusp
compute_θmin(model::AbstractInterpolatorProfile) = exp(first(first(model.itp.ranges)))
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
struct LogInterpolatorProfile{T, P <: AbstractProfile{T}, I1, C} <: AbstractInterpolatorProfile{T}
    model::P
    itp::I1
    cosmo::C
end


function LogInterpolatorProfile(model::AbstractProfile, itp)
    return LogInterpolatorProfile(model, itp, model.cosmo)  # use wrapped cosmology
end

# forward the interpolator calls to the wrapped interpolator
# IMPORTANT: for backwards compat, interpolator internal order is θ, z, mass
# which DIFFERS from the rest of the code which is (θ, mass, z, α, δ)
# should fix this at some point
@inline (ip::LogInterpolatorProfile)(θ, Mh_Msun, z) = exp(ip.itp(log(θ), z, log10(Mh_Msun)))

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
                        workspace::CarClenshawCurtisProfileWorkspace, model, Mh, z, α₀, δ₀, 
                        θmax, normalization=1) where T

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
                             T(normalization * model(θ, Mh, z)),
                             zero(T))
        end
    end
end

# fall back to generic profile painter if no specialized painter is defined for the model
function profile_paint!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, 
                        workspace::CarClenshawCurtisProfileWorkspace, model, 
                        Mh, z, α₀, δ₀, θmax, normalization=1) where T
    profile_paint_generic!(m, workspace, model, Mh, z, α₀, δ₀, θmax, normalization)
end


function profile_paint_generic!(m::Enmap{T, 2, Matrix{T}, Gnomonic{T}}, 
                                workspace::GnomonicProfileWorkspace, model, 
                                Mh, z, α₀, δ₀, θmax, normalization=1) where T

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
                             normalization * model(θ, Mh, z),
                             zero(T))
        end
    end
end

# fall back to generic profile painter if no specialized painter is defined for the model
function profile_paint!(m::Enmap{T, 2, Matrix{T}, Gnomonic{T}}, 
                        workspace::GnomonicProfileWorkspace, model, Mh, z, α₀, δ₀, 
                        θmax, normalization=1) where T
    profile_paint_generic!(m, model, workspace, Mh, z, α₀, δ₀, θmax, normalization)
end


function profile_paint_generic!(m::HealpixMap{T, RingOrder}, w::HealpixSerialProfileWorkspace, 
        model, Mh, z, α₀, δ₀, θmax, normalization=1) where T
    ϕ₀ = α₀
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    θmin = compute_θmin(model)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(clamp(1 - d² / 2, -one(T), one(T)))
        θ = max(θmin, θ)  # clamp to minimum θ
        m.pixels[ir] += ifelse(θ < θmax, 
                                    normalization * model(θ, Mh, z),
                                    zero(T))
    end
end

# fall back to generic profile painter if no specialized painter is defined for the model
function profile_paint!(m::HealpixMap{T, RingOrder}, w::HealpixSerialProfileWorkspace, model, 
                        Mh, z, α₀, δ₀, θmax, normalization=1) where T
    profile_paint_generic!(m, w, model, Mh, z, α₀, δ₀, θmax, normalization)
end


# paint the the sources in the given range of indices
function paintrange!(irange::AbstractUnitRange, m, workspace, model, masses, redshifts, αs, δs)
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        Mh = masses[i]
        z = redshifts[i]
        θmax_ = compute_θmax(model, Mh * XGPaint.M_sun, z)
        profile_paint!(m, workspace, model, Mh, z, α₀, δ₀, θmax_)
    end
end


_fillzero!(m) = fill!(m, zero(eltype(m)))
_fillzero!(m::HealpixMap) = fill!(m.pixels, zero(eltype(m)))

# paint! is threaded by default
function paint!(m, workspace, model, masses, redshifts, αs, δs; 
                zerobeforepainting=true)
    
    zerobeforepainting && _fillzero!(m)

    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize)

    if N_sources < 2Threads.nthreads()  # don't thread if there are not many sources
        return paintrange!(1:N_sources, m, wrapserialworkspace(workspace, 1), 
            model, masses, redshifts, αs, δs)
    end

    Threads.@threads for ti in 1:Threads.nthreads()
        chunk_i = 2ti
        i1, i2 = chunks[chunk_i]
        paintrange!(i1:i2, m, wrapserialworkspace(workspace, ti), 
            model, masses, redshifts, αs, δs)
    end

    Threads.@threads for ti in 1:Threads.nthreads()
        chunk_i = 2ti - 1
        i1, i2 = chunks[chunk_i]
        paintrange!(i1:i2, m, wrapserialworkspace(workspace, ti), 
            model, masses, redshifts, αs, δs)
    end
end


# for kSZ, we need to extend paintrange! and paint! to take in a velocity

# paint the the sources in the given range
function paintrange!(irange::AbstractUnitRange, m, workspace, model, 
                     masses, redshifts, αs, δs, proj_v_over_c)
    for i in irange
        θmax = compute_θmax(model, masses[i] * XGPaint.M_sun, redshifts[i])
        profile_paint!(m, workspace, model, 
            masses[i], redshifts[i], αs[i], δs[i], θmax, proj_v_over_c[i])
    end
end


# extend general paint! to take in a projected velocity
function paint!(m, workspace, model, masses, redshifts, αs, δs, proj_v_over_c; 
        zerobeforepainting=true)
    zerobeforepainting && _fillzero!(m)

    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize)

    if N_sources < 2Threads.nthreads()  # don't thread if there are not many sources
        return paintrange!(1:N_sources, m, wrapserialworkspace(workspace, 1), 
            model, masses, redshifts, αs, δs, proj_v_over_c)
    end

    Threads.@threads for ti in 1:Threads.nthreads()
        chunk_i = 2ti
        i1, i2 = chunks[chunk_i]
        paintrange!(i1:i2, m, wrapserialworkspace(workspace, ti), 
            model, masses, redshifts, αs, δs, proj_v_over_c)
    end

    Threads.@threads for ti in 1:Threads.nthreads()
        chunk_i = 2ti - 1
        i1, i2 = chunks[chunk_i]
        paintrange!(i1:i2, m, wrapserialworkspace(workspace, ti), 
            model, masses, redshifts, αs, δs, proj_v_over_c)
    end
end


# # serial version of the paint function, mostly for debugging
# function paint!(m, workspace::HealpixSerialProfileWorkspace, model, masses, redshifts, αs, δs; 
#         zerobeforepainting=true)
#     zerobeforepainting && _fillzero!(m)
#     return paintrange!(1:length(masses), m, wrapserialworkspace(workspace, 1), 
#         model, masses, redshifts, αs, δs)
# end
