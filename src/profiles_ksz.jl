
import PhysicalConstants.CODATA2018 as constants
using Unitful
const M_sun = 1.98847e30u"kg"
const P_e_factor = constants.σ_e / (constants.m_e * constants.c_0^2)
const T_cmb =  2.725 * u"K"
using Cosmology
using QuadGK
using DelimitedFiles
using Interpolations
using NPZ


table = npzread("/fs/lustre/cita/zack/projects/SZpack.v1.1.1/szpack_interp_4d.npz")
nu_vector = LinRange(log(5.680062373019096*1e9),log(852.0093559528645*1e9),500)
temp_vector = LinRange(1.0e-3,75.0,100)
vel_vector = LinRange(1.0e-9, 1.0e-1,100)
mu_vector = LinRange(-1,1,50)
szpack_interp_ksz = scale(Interpolations.interpolate(table, BSpline(Cubic(Line(OnGrid())))), (temp_vector), (vel_vector), (mu_vector), (nu_vector))

table_T0 = npzread("/fs/lustre/cita/zack/projects/SZpack.v1.1.1/szpack_interp_4d_0.npz")
szpack_interp_T0 = scale(Interpolations.interpolate(table_T0[1,:,:,:], BSpline(Cubic(Line(OnGrid())))), (vel_vector), (mu_vector), (nu_vector))


struct Battaglia16KinematicSZProfile{T,C} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
    X::T  # X = 2.6408 corresponding to frequency 150 GHz
end

function Battaglia16KinematicSZProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, x::T=2.6408) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    X = x
    return Battaglia16KinematicSZProfile(f_b, cosmo, X)
end


function dimensionless_P_profile_los_ksz(𝕡::Battaglia16KinematicSZProfile{T}, M_200, z, r) where T
    par = get_params(𝕡, M_200, z)
    R_200 = R_Δ(𝕡, M_200, z, 200)
    x = r / angular_size(𝕡, R_200, z)
    return par.P₀ * _tsz_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
end


"""Line-of-sight integrated electron pressure"""
P_e_los_ksz(𝕡, M_200, z, r) = 0.5176 * P_th_los_ksz(𝕡, M_200, z, r)

"""Line-of-sight integrated thermal pressure"""
P_th_los_ksz(𝕡, M_200, z, r) = constants.G * M_200 * 200 * ρ_crit(𝕡, z) * 
    𝕡.f_b / 2 * dimensionless_P_profile_los_ksz(𝕡, M_200, z, r)


function compton_y_ksz(𝕡, M_200, z, r)
    return P_e_los_ksz(𝕡, M_200, z, r) * P_e_factor
end


function SZpack_ksz(𝕡, M_200, z, r; vel=3e3, τ=0.01, mu = 1.0, showT=true)
    """
    Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
    """
    X = 𝕡.X
    T_e = T_vir_calc(𝕡, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(X_to_nu(X)))
    vel = abs(ustrip(vel/uconvert(u"km/s",constants.c_0))) # need to take absolute value of v to make sure v is within bounds of interpolator
        
    # Term 1
    dI_1 = uconvert(u"kg*s^-2",(szpack_interp_ksz(t, vel, mu, nu)*u"MJy/sr" - szpack_interp_T0(vel, mu, nu)*u"MJy/sr")/τ)
    y = XGPaint.compton_y_ksz(𝕡, M_200, z, r)
    I_1 = uconvert(u"kg*s^-2",y * (dI_1/(θ_e)))
    
    # Term 2
    dI_2 = uconvert(u"kg*s^-2",(szpack_interp_T0(vel, mu, nu)*u"MJy/sr")/τ)
    tau = 0.01 #XGPaint.tau_ksz(𝕡, M_200, z, r)
    I_2 = uconvert(u"kg*s^-2",dI_2 * tau)
    
    I = I_1 + I_2
    T = I/uconvert(u"kg*s^-2",abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X)/(constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))

    if showT==true
        return abs(T)
    else
        return I
    end
end


function profile_grid_ksz(𝕡::AbstractGNFW{T}; N_z=256, N_logM=256, N_logθ=512, z_min=1e-3, z_max=5.0, 
              logM_min=11, logM_max=15.7, logθ_min=-16.5, logθ_max=2.5, N_v=256, v_min=-1353.6917, v_max=1276.7216) where T #logM_max=15.7

    logθs = LinRange(logθ_min, logθ_max, N_logθ)
    redshifts = LinRange(z_min, z_max, N_z)
    logMs = LinRange(logM_min, logM_max, N_logM)
    velocities = LinRange(v_min, v_max, N_v)

    return profile_grid_ksz(𝕡, logθs, redshifts, logMs, velocities)
end


function profile_grid_ksz(𝕡::AbstractGNFW{T}, logθs, redshifts, logMs, velocities) where T

    N_logθ, N_z, N_logM, N_vels = length(logθs), length(redshifts), length(logMs), length(velocities)
    A = zeros(T, (N_logθ, N_z, N_logM, N_vels))

    Threads.@threads for im in 1:N_logM
        println("Completed Halo Mass $im")
        logM = logMs[im]
        M = 10^(logM) * M_sun
        for (iz, z) in enumerate(redshifts)
            for iθ in 1:N_logθ
                θ = exp(logθs[iθ])
                for (iv, vel) in enumerate(velocities)
                    szp = SZpack_ksz(𝕡, M, z, θ; vel=vel)
                    A[iθ, iz, im, iv] = max(zero(T), szp)
                end
            end
        end
    end

    return logθs, redshifts, logMs, velocities, A
end


function profile_paint_ksz!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, p,
                        α₀, δ₀, psa::CarClenshawCurtisProfileWorkspace, 
                        sitp, z, Ms, vel, θmax) where T
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
                                 #sign * exp(sitp(log(θ), z, log10(Ms))),
                                 exp(sitp(log(θ), z, log10(Ms), vel)),
                                   zero(T))
        end
    end
end


"""Helper function to build a tSZ interpolator"""
function build_interpolator_ksz(model::AbstractGNFW; cache_file::String="",
                            N_logθ=512, pad=256, overwrite=true, verbose=true)

    if overwrite || (isfile(cache_file) == false)
        verbose && print("Building new interpolator from model.\n")
        rft = RadialFourierTransform(n=N_logθ, pad=pad)
        logθ_min, logθ_max = log(minimum(rft.r)), log(maximum(rft.r))
        prof_logθs, prof_redshift, prof_logMs, prof_velocity, prof_y = profile_grid_ksz(model; 
            N_logθ=N_logθ, logθ_min=logθ_min, logθ_max=logθ_max)
        if length(cache_file) > 0
            verbose && print("Saving new interpolator to $(cache_file).\n")
            save(cache_file, Dict("prof_logθs"=>prof_logθs, 
                "prof_redshift"=>prof_redshift, "prof_logMs"=>prof_logMs, "prof_velocity"=>prof_velocity, "prof_y"=>prof_y))
        end
    else
        print("Found cached Battaglia profile model. Loading from disk.\n")
        model_grid = load(cache_file)
        prof_logθs, prof_redshift, prof_logMs, prof_velocity, prof_y = model_grid["prof_logθs"], 
            model_grid["prof_redshift"], model_grid["prof_logMs"], model_grid["prof_velocity"], model_grid["prof_y"]
    end
    
    itp = Interpolations.interpolate(log.(prof_y), BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, prof_logθs, prof_redshift, prof_logMs, prof_velocity)
    return sitp
end


function profile_paint_ksz!(m::HealpixMap{T, RingOrder}, p,
            α₀, δ₀, w::HealpixProfileWorkspace, z, Mh, vel, θmax) where T
    ϕ₀ = α₀
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    sitp = w.profile_real_interp
    
    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(1 - d² / 2)
        θ = max(w.θmin, θ)  # clamp to minimum θ
        m.pixels[ir] += ifelse(θ < θmax, 
                                   #sign * exp(sitp(log(θ), z, log10(Mh))),
                                   exp(sitp(log(θ), z, log10(Mh), vel)),
                                    zero(T))
    end
end


# for rectangular pixelizations

# multi-halo painting utilities
function paint_ksz!(m, p::XGPaint.AbstractProfile, psa, sitp, 
                masses::AV, redshifts::AV, αs::AV, δs::AV, velocities::AV, irange::AbstractUnitRange) where AV
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        v = velocities[i]
        θmax_ = θmax(p, mh * XGPaint.M_sun, z)
        profile_paint_ksz!(m, p, α₀, δ₀, psa, sitp, z, mh, v, θmax_)
    end
end


function paint_ksz!(m, p::XGPaint.AbstractProfile, psa, sitp, masses::AV, 
                        redshifts::AV, αs::AV, δs::AV, velocities::AV)  where AV
    fill!(m, 0)
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_ksz!(m, p, psa, sitp, masses, redshifts, αs, δs,velocities, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_ksz!(m, p, psa, sitp, masses, redshifts, αs, δs, velocities, i1:i2)
    end
end

