# ──────────────────────────────────────────────────────────────────────────────
#  FAST 1‑D SPLINE FOR ARNAUD ET AL. (2010)  (θ , M₅₀₀ , z) → y
# ──────────────────────────────────────────────────────────────────────────────

@inline function _a10_build_table(model::Arnauld10ThermalSZProfile{T};
                                  Nx::Int = 256,
                                  x_min   = 1e-12,
                                  x_max   = 1e5) where T
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
