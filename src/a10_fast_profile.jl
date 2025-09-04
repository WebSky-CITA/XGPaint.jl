# ──────────────────────────────────────────────────────────────────────────────
#  FAST 1‑D SPLINE FOR ARNAUD ET AL. (2010)  (θ , M₅₀₀ , z) → y
# ──────────────────────────────────────────────────────────────────────────────

@inline function _a10_build_table(model::Arnauld10ThermalSZProfile{T};
                                  Nx::Int = 2048,
                                  x_min   = 1e-12,
                                  x_max   = 1e12) where T
    logx = range(log(x_min), log(x_max); length = Nx)   # StepRangeLen
    x    = exp.(logx)

    par  = get_params(model, 1e14 * M_sun, 0.0)
    F    = par.P₀ .* _nfw_profile_los_quadrature.(x,
                                                  par.xc, par.α, par.β, par.γ)

    return Interpolations.scale(
               Interpolations.interpolate(log.(F),
                                          BSpline(Cubic(Line(OnGrid())))),
               logx)                                   # now valid
end


struct A10ThetaProfile{T,I,C} <: AbstractInterpolatorProfile{T}
    itp  :: I                                  # spline in log x
    base :: Arnauld10ThermalSZProfile{T,C}     # keeps B & cosmology
end

A10ThetaProfile(model::Arnauld10ThermalSZProfile{T}; kwargs...) where T =
    A10ThetaProfile(_a10_build_table(model; kwargs...), model)


# ── evaluation ───────────────────────────────────────────────────────────────
@inline function (p::A10ThetaProfile)(θ::Real, M500::Real, z::Real)
    m    = p.base
    R500 = R_Δ(m, M500 * M_sun, z, 500) / m.B^(1/3)
    θ500 = angular_size(m, R500, z)
    x    = θ / θ500

    Fval  = exp(p.itp(log(x)))
    Pe500 = A10_normalization(m, M500 * M_sun, z; B = m.B)

    return Pe500 * R500 * Fval * P_e_factor + 0   # +0 strips Unitful units
end


# ── helpers so painters behave exactly as before ─────────────────────────────
compute_θmin(p::A10ThetaProfile) = exp(first(first(p.itp.ranges)))

# forward R_Δ and angular_size (and anything else) to the wrapped base model
R_Δ(p::A10ThetaProfile, args...)        = R_Δ(p.base,        args...)
angular_size(p::A10ThetaProfile, args...) = angular_size(p.base, args...)

# compute_θmax used in paintrange!/paint!
function compute_θmax(p::A10ThetaProfile{T}, MΔ, z; kwargs...) where T
    compute_θmax(p.base, MΔ, z; kwargs...)   # delegate
end

# make .cosmo and .B visible if accessed directly
Base.getproperty(p::A10ThetaProfile, s::Symbol) =
    (s === :cosmo || s === :B) ? getfield(p.base, s) : getfield(p, s)

Base.show(io::IO, p::A10ThetaProfile) =
    print(io, "A10ThetaProfile( grid = ", length(first(p.itp.ranges)), " pts )")


# ── overload build_interpolator just for A10, leaving others untouched ───────
function build_interpolator(model::Arnauld10ThermalSZProfile; kwargs...)
    A10ThetaProfile(model; kwargs...)
end


############ 2D grid in (θ, θ500) + beam on θ (FHT) + cubic–cubic 2D spline ############
# Assumes:
# - A10ThetaProfile is defined (your code above)
# - RadialFourierTransform, real2harm, harm2real, gaussbeam are available
# - Interpolations.jl is loaded

# 1) Build the raw 2D shape grid: F(θ, θ500) := F(θ/θ500) on uniform log axes
# ---------------------------------------------------------------------
# 2D grid in (logθ, logθ500) for the universal A10 shape F(u=θ/θ500)
# ---------------------------------------------------------------------

function profile_grid_θ_θ500(shape::A10ThetaProfile{T};
                             N_logθ::Int=256, N_logθ500::Int=128,
                             logθ_min::Real=-16.5,  logθ_max::Real=2.5,
                             logθ500_min::Real=-8.0, logθ500_max::Real=-2.0) where T
    logθs     = LinRange(logθ_min,   logθ_max,   N_logθ)
    logθ500s  = LinRange(logθ500_min, logθ500_max, N_logθ500)
    return profile_grid_θ_θ500(shape, logθs, logθ500s)
end

function profile_grid_θ_θ500(shape::A10ThetaProfile{T},
                             logθs::AbstractVector,
                             logθ500s::AbstractVector) where T
    N_logθ     = length(logθs)
    N_logθ500  = length(logθ500s)
    A          = zeros(T, (N_logθ, N_logθ500))

    umin = compute_θmin(shape)  # minimal u supported by the 1D spline

    Threads.@threads :static for j in 1:N_logθ500
        θ500 = exp(logθ500s[j])
        for iθ in 1:N_logθ
            θ = exp(logθs[iθ])
            u = max(umin, θ/θ500)
            # mirror your original 'max(zero(T), ...)' safety
            A[iθ, j] = max(zero(T), T(exp(shape.itp(log(u)))))
        end
    end

    return logθs, logθ500s, A
end

# 2) Beam-convolve along θ (first axis) using the same padded FHT pipeline as your old code
function transform_profile_grid_θ!(y_prof_grid, rft, lbeam)
    pad_val      = rft.pad
    Nθ, Nθ500    = size(y_prof_grid)

    padded = zeros(eltype(y_prof_grid), Nθ + 2*pad_val, Nθ500)
    padded[pad_val+1 : pad_val+Nθ, :] .= y_prof_grid

    for j in axes(padded, 2)
        rprof = padded[:, j]             # COPY (no views), like original
        lprof = real2harm(rft, rprof)
        lprof .*= lbeam
        reverse!(lprof)                  # same convention as original
        rprof_t = harm2real(rft, lprof)
        padded[:, j] .= rprof_t
    end

    y_prof_grid .= padded[pad_val+1 : pad_val+Nθ, :]
    return y_prof_grid
end


# 3) Clean up small negatives after numerics (same logic, now 2D: θ as first axis)
"prune a 2D (θ × θ500) grid for negatives, extrapolate forward along θ"
function cleanup_negatives_θ!(A::AbstractMatrix)
    for j in axes(A, 2)                   # θ500 index
        extrapolating = false
        fact = 1.0
        for k in axes(A, 1)               # θ index
            if A[k, j] <= 0
                extrapolating = true
                fact = A[k-1, j] / A[k-2, j]
            end
            if extrapolating
                A[k, j] = max(fact * A[k-1, j], nextfloat(0.0))
            end
        end
    end
    return A
end

# === 2D beamed A10 wrapper (type must be defined before methods) ===
struct BeamConvolvedA10_2D{I,B}
    itp2d :: I   # FilledExtrapolation( ScaledInterpolation( ... ) )
    base  :: B   # Arnauld10ThermalSZProfile (geometry & normalization)
end

# === Cached drop-in model: same call (θ, M, z); avoids per-pixel cosmology ===
mutable struct CachedBeamConvolvedA10_2D{I,B,T}
    itp2d   :: I
    base    :: B
    last_M  :: T
    last_z  :: T
    logθ500 :: Float64
    Afac    :: T
    θmin    :: T
    filled  :: Bool
end

function CachedBeamConvolvedA10_2D(p::BeamConvolvedA10_2D{I,B}) where {I,B}
    T = Float64
    θmin = compute_θmin(p)
    return CachedBeamConvolvedA10_2D{I,B,T}(p.itp2d, p.base, NaN, NaN, NaN, zero(T), θmin, false)
end

# Delegate helpers so existing code keeps working
R_Δ(p::CachedBeamConvolvedA10_2D, args...)          = R_Δ(p.base, args...)
angular_size(p::CachedBeamConvolvedA10_2D, args...) = angular_size(p.base, args...)
Base.getproperty(p::CachedBeamConvolvedA10_2D, s::Symbol) =
    (s === :cosmo || s === :B) ? getfield(p.base, s) : getfield(p, s)

compute_θmin(p::CachedBeamConvolvedA10_2D) = p.θmin
function compute_θmax(p::CachedBeamConvolvedA10_2D, M_Δ, z; FWHM=10, mult=4)
    r      = R_Δ(p, M_Δ, z)
    theta1 = 2 * FWHM * (π/10800)
    theta2 = mult * angular_size(p, r, z)
    return oftype(theta2, max(theta1, theta2))
end

@inline function (p::CachedBeamConvolvedA10_2D{I,B,T})(θ::Real, M500::Real, z::Real) where {I,B,T}
    if !p.filled || M500 != p.last_M || z != p.last_z
        m     = p.base
        R500  = R_Δ(m, M500 * M_sun, z, 500) / m.B^(1/3)
        θ500  = angular_size(m, R500, z)
        p.logθ500 = log(Float64(θ500))
        p.Afac    = A10_normalization(m, M500 * M_sun, z; B = m.B) * R500 * P_e_factor
        p.last_M  = T(M500)
        p.last_z  = T(z)
        p.filled  = true
    end
    θeff = max(p.θmin, T(θ))
    return p.Afac * p.itp2d(log(Float64(θeff)), p.logθ500) + 0
end


# 1) Make the wrapper callable: model(θ, M500, z)
@inline function (p::BeamConvolvedA10_2D)(θ::Real, M500::Real, z::Real)
    m     = p.base
    R500  = R_Δ(m, M500 * M_sun, z, 500) / m.B^(1/3)
    θ500  = angular_size(m, R500, z)
    Afac  = A10_normalization(m, M500 * M_sun, z; B = m.B) * R500 * P_e_factor
    return Afac * p.itp2d(log(Float64(θ)), log(Float64(θ500))) + 0
end

# 2) Geometry forwarders (so compute_θmax etc. can delegate)
R_Δ(p::BeamConvolvedA10_2D, args...)          = R_Δ(p.base, args...)
angular_size(p::BeamConvolvedA10_2D, args...) = angular_size(p.base, args...)

# 3) Expose cosmology/boost like the 1-D wrapper (fixes ρ_crit → .cosmo access)
Base.getproperty(p::BeamConvolvedA10_2D, s::Symbol) =
    (s === :cosmo || s === :B) ? getfield(p.base, s) : getfield(p, s)


# ──────────────────────────────────────────────────────────────────────────────

# helper to reach ranges through FilledExtrapolation/ScaledInterpolation stacks
@inline _ranges2d(itp) = hasproperty(itp, :ranges) ? itp.ranges : _ranges2d(getfield(itp, :itp))

# θ_min from the first (logθ) axis of the 2-D interpolator
function compute_θmin(p::BeamConvolvedA10_2D)
    logθ_axis = first(_ranges2d(p.itp2d))
    return exp(first(logθ_axis))
end

# θ_max: same rule as original, but forward geometry to the wrapped base model
function compute_θmax(p::BeamConvolvedA10_2D, M_Δ, z; FWHM=10, mult=4)
    r      = R_Δ(p, M_Δ, z)                       # forwards to p.base
    theta1 = 2 * FWHM * (π/10800)                 # 2×FWHM in radians
    theta2 = mult * angular_size(p, r, z)         # forwards to p.base
    return oftype(theta2, max(theta1, theta2))    # preserve element type
end

# 2D variant: (logθ, logθ500, A[θ, θ500]) → interpolator of logR vs logθ500
function build_max_paint_logradius(logθs::AbstractVector,
                                   logθ500s::AbstractVector,
                                   A::AbstractArray{T,2};
                                   rtol::Real = 1e-2) where {T}

    N_logθ     = length(logθs)
    N_logθ500  = length(logθ500s)
    logRs      = zeros(T, N_logθ500)     # output: log R(θ500)
    dF_r       = zeros(T, N_logθ)        # cumulative integral along θ

    for j in 1:N_logθ500
        s = zero(T)
        # integrate F(θ, θ500) * θ dθ (trapezoid) along θ
        for iθ in 1:(N_logθ-1)
            θ₁ = exp(logθs[iθ])
            θ₂ = exp(logθs[iθ+1])
            f₁ = A[iθ,   j] * θ₁
            f₂ = A[iθ+1, j] * θ₂
            s += (θ₂ - θ₁) * (f₁ + f₂) / 2
            dF_r[iθ] = s
        end

        if s <= 0
            # degenerate column; fall back to smallest θ bin
            logRs[j] = logθs[1]
            continue
        end

        threshold = (one(T) - T(rtol)) * s
        # find last θ where cumulative < threshold, like original
        for iθ in (N_logθ-1):-1:1
            if dF_r[iθ] < threshold
                logRs[j] = min(logθs[iθ], log(π))
                break
            end
        end
    end

    # return a 1D cubic spline of logR vs logθ500 (same basis as original)
    itp = Interpolations.interpolate(
        logRs,
        Interpolations.BSpline(Interpolations.Cubic(Interpolations.Line(Interpolations.OnGrid())))
    )
    return Interpolations.scale(itp, logθ500s)
end
