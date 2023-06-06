
import PhysicalConstants.CODATA2018 as constants
using Unitful
const M_sun = 1.98847e30u"kg"
const P_e_factor = constants.σ_e / (constants.m_e * constants.c_0^2)
using Cosmology
using QuadGK



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

struct BattagliaProfile{T,C} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
end

struct BreakModel{T,C} <: AbstractGNFW{T}
    f_b::T
    cosmo::C
    alpha_break::T
    M_break::T
end

function BattagliaProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    return BattagliaProfile(f_b, cosmo)
end

abstract type AbstractPaintingProblem{T} end


function BreakModel(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, alpha_break::T=1.5, M_break::T=2.0*10^14) where {T <: Real}
    #alpha_break = 1.486 from Shivam P paper by Nate's sleuthing
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    return BreakModel(f_b, cosmo, alpha_break, M_break)
end

const ρ_crit_factor = uconvert(u"kg/m^3", 3u"km^2*Mpc^-2*s^-2" / (8π * constants.G))


function ρ_crit(𝕡, z)
    H_z = H(𝕡.cosmo, z)
    return uconvert(u"kg/m^3", 3H_z^2 / (8π * constants.G))
end

function R_Δ(𝕡, M_Δ, z, Δ=200)
    return ∛(M_Δ / (4π/3 * Δ * ρ_crit(𝕡, z)))
end

function angular_size(𝕡::AbstractProfile{T}, physical_size, z) where T
    d_A = angular_diameter_dist(𝕡.cosmo, z)

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
function _tsz_profile_los_quadrature(x, xc, α, β, γ; zmax=1e5, rtol=eps(), order=9)
    x² = x^2
    scale = 1e9
    integral, err = quadgk(y -> scale * generalized_nfw(√(y^2 + x²), xc, α, β, γ),
                      0.0, zmax, rtol=rtol, order=order)
    return 2integral / scale
end

function dimensionless_P_profile_los(𝕡::BattagliaProfile{T}, M_200, z, r) where T
    par = get_params(𝕡, M_200, z)
    R_200 = R_Δ(𝕡, M_200, z, 200)
    x = r / angular_size(𝕡, R_200, z)
    return par.P₀ * _tsz_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
end

function dimensionless_P_profile_los(𝕡::BreakModel{T}, M_200, z, r) where T
    par = get_params(𝕡, M_200, z)
    R_200 = R_Δ(𝕡, M_200, z, 200)
    x = r / angular_size(𝕡, R_200, z)
    if M_200 < 𝕡.M_break * M_sun
        return par.P₀ * (M_200/(𝕡.M_break*M_sun))^𝕡.alpha_break * _tsz_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
    else
        return par.P₀ * _tsz_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
    end
end

"""Line-of-sight integrated electron pressure"""
P_e_los(𝕡, M_200, z, r) = 0.5176 * P_th_los(𝕡, M_200, z, r)

"""Line-of-sight integrated thermal pressure"""
P_th_los(𝕡, M_200, z, r) = constants.G * M_200 * 200 * ρ_crit(𝕡, z) * 
    𝕡.f_b / 2 * dimensionless_P_profile_los(𝕡, M_200, z, r)

function compton_y(𝕡, M_200, z, r)
    return P_e_los(𝕡, M_200, z, r) * P_e_factor
end

function profile_grid(𝕡::AbstractGNFW{T}; N_z=256, N_logM=256, N_logθ=512, z_min=1e-3, z_max=5.0, 
              logM_min=11, logM_max=15.7, logθ_min=-16.5, logθ_max=2.5) where T

    logθs = LinRange(logθ_min, logθ_max, N_logθ)
    redshifts = LinRange(z_min, z_max, N_z)
    logMs = LinRange(logM_min, logM_max, N_logM)

    return profile_grid(𝕡, logθs, redshifts, logMs)
end

function profile_grid(𝕡::AbstractGNFW{T}, logθs, redshifts, logMs) where T

    N_logθ, N_z, N_logM = length(logθs), length(redshifts), length(logMs)
    A = zeros(T, (N_logθ, N_z, N_logM))

    Threads.@threads for im in 1:N_logM
        logM = logMs[im]
        M = 10^(logM) * M_sun
        for (iz, z) in enumerate(redshifts)
            for iθ in 1:N_logθ
                θ = exp(logθs[iθ])
                y = compton_y(𝕡, M, z, θ)
                A[iθ, iz, im] = max(zero(T), y)
            end
        end
    end

    return logθs, redshifts, logMs, A
end


# get angular size in radians of radius to stop at
function θmax(𝕡::AbstractProfile{T}, M_Δ, z; mult=4) where T
    r = R_Δ(𝕡, M_Δ, z)
    return T(mult * angular_size(𝕡, r, z))
end

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


function profile_paint!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, 
                        α₀, δ₀, psa::CarClenshawCurtisProfileWorkspace, 
                        sitp, z, Ms, θmax) where T

    # get indices of the region to work on
    i1, j1 = sky2pix(m, α₀ - θmax, δ₀ - θmax)
    i2, j2 = sky2pix(m, α₀ + θmax, δ₀ + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))

    x₀ = cos(δ₀) * cos(α₀)
    y₀ = cos(δ₀) * sin(α₀) 
    z₀ = sin(δ₀)

    @inbounds for j in j_start:j_stop
        for i in i_start:i_stop
            x₁ = psa.cos_δ[i,j] * psa.cos_α[i,j]
            y₁ = psa.cos_δ[i,j] * psa.sin_α[i,j]
            z₁ = psa.sin_δ[i,j]
            d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
            θ =  acos(1 - d² / 2)
            m[i,j] += ifelse(θ < θmax, 
                             exp(sitp(log(θ), z, log10(Ms))),
                             zero(T))
        end
    end
end


function profile_paint!(m::Enmap{T, 2, Matrix{T}, Gnomonic{T}}, 
            α₀, δ₀, psa::GnomonicProfileWorkspace, sitp, z, Ms, θmax) where T

    # get indices of the region to work on
    i1, j1 = sky2pix(m, α₀ - θmax, δ₀ - θmax)
    i2, j2 = sky2pix(m, α₀ + θmax, δ₀ + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))

    x₀ = cos(δ₀) * cos(α₀)
    y₀ = cos(δ₀) * sin(α₀) 
    z₀ = sin(δ₀)

    @inbounds for j in j_start:j_stop
        for i in i_start:i_stop
            x₁ = psa.cos_δ[i,j] * psa.cos_α[i,j]
            y₁ = psa.cos_δ[i,j] * psa.sin_α[i,j]
            z₁ = psa.sin_δ[i,j]
            d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
            θ =  acos(1 - d² / 2)
            m[i,j] += ifelse(θ < θmax, 
                             exp(sitp(log(θ), z, log10(Ms))),
                             zero(T))
        end
    end
end


function profile_paint!(m::HealpixMap{T, RingOrder}, 
            α₀, δ₀, w::HealpixProfileWorkspace, z, Mh, θmax) where T
    ϕ₀ = α₀
    θ₀ = π/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    sitp = w.profile_real_interp
    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(1 - d² / 2)
        θ = max(w.θmin, θ)  # clamp to minimum θ
        m.pixels[ir] += ifelse(θ < θmax, 
                                    exp(sitp(log(θ), z, log10(Mh))),
                                    zero(T))
    end
end





# for rectangular pixelizations

# multi-halo painting utilities
function paint!(m, p::XGPaint.AbstractProfile, psa, sitp, 
                masses::AV, redshifts::AV, αs::AV, δs::AV, irange::AbstractUnitRange) where AV
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax_ = θmax(p, mh * XGPaint.M_sun, z)
        profile_paint!(m, α₀, δ₀, psa, sitp, z, mh, θmax_)
    end
end

function paint!(m, p::XGPaint.AbstractProfile, psa, sitp, masses::AV, 
                        redshifts::AV, αs::AV, δs::AV)  where AV
    m .= 0.0
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2)
    end
end



