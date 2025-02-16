

function read_szpack_table(filename)
    table = readdlm(filename)
    nu_vector = LinRange(log(5.680062373019096*1e9),log(852.0093559528645*1e9),3000)
    temp_vector = LinRange(1.0e-3,75.0,100)
    szpack_interp = scale(Interpolations.interpolate(table, BSpline(Cubic(Line(OnGrid())))), (temp_vector), (nu_vector))
    return szpack_interp
end

function X_to_nu(X)
    return (X*constants.k_B*T_cmb)/constants.h
end

nu_to_X(nu) = (constants.h*nu)/(constants.k_B*T_cmb)

struct Battaglia16SZPackProfile{T,C,TSZ, ITP1, ITP2} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
    X::T  # X = 2.6408 corresponding to frequency 150 GHz
    𝕡_tsz::TSZ
    y_interp::ITP1
    szpack_interp::ITP2
    τ::T
end

function Battaglia16SZPackProfile(𝕡_tsz, y_interp, x::T, τ=0.01; Omega_c=0.2589, 
        Omega_b=0.0486, h=0.6774, table_filename=rsz_szpack_table_filename()) where T
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    X = x
    szpack_interp = read_szpack_table(table_filename)
    return Battaglia16SZPackProfile(f_b, cosmo, X, 𝕡_tsz, y_interp, szpack_interp, τ)
end

function SZpack(𝕡, M_200, z, r; τ=0.01, showT=true)
    """
    Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
    """
    X = 𝕡.X
    T_e = T_vir_calc(𝕡, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)

    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(X_to_nu(X)))
    dI = uconvert(u"kg*s^-2",𝕡.szpack_interp(t, nu)*u"MJy/sr")
    
    y = XGPaint.compton_y(𝕡.𝕡_tsz, M_200, z, r)
    I = uconvert(u"kg*s^-2",y * (dI/(τ * θ_e)))
    T = I/uconvert(u"kg*s^-2",abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X)/(constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))

    if showT==true
        return abs(T)
    else
        return I
    end
end


function profile_grid_szp(𝕡::AbstractGNFW{T}; N_z=256, N_logM=256, N_logθ=512, z_min=1e-3, z_max=5.0, 
              logM_min=11, logM_max=15.7, logθ_min=-16.5, logθ_max=2.5) where T

    logθs = LinRange(logθ_min, logθ_max, N_logθ)
    redshifts = LinRange(z_min, z_max, N_z)
    logMs = LinRange(logM_min, logM_max, N_logM)

    return profile_grid_szp(𝕡, logθs, redshifts, logMs)
end


function profile_grid_szp(𝕡::AbstractGNFW{T}, logθs, redshifts, logMs) where T

    N_logθ, N_z, N_logM = length(logθs), length(redshifts), length(logMs)
    A = zeros(T, (N_logθ, N_z, N_logM))

    Threads.@threads for im in 1:N_logM
        logM = logMs[im]
        M = 10^(logM) * M_sun
        for (iz, z) in enumerate(redshifts)
            for iθ in 1:N_logθ
                θ = exp(logθs[iθ])
                szp = SZpack(𝕡, M, z, θ)
                A[iθ, iz, im] = max(zero(T), szp)
            end
        end
    end

    return logθs, redshifts, logMs, A
end


function I_to_T_mult_factor(X)
    return 1/(abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X) / 
        (constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))
end

function profile_paint_szp!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, 
                        p::Battaglia16SZPackProfile, 
                        α₀, δ₀, psa::CarClenshawCurtisProfileWorkspace, 
                        z, Mh, θmax) where T
    # get indices of the region to work on
    i1, j1 = sky2pix(m, α₀ - θmax, δ₀ - θmax)
    i2, j2 = sky2pix(m, α₀ + θmax, δ₀ + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))
    θmin = exp(first(first(interp_model.itp.ranges)))

    # needs mass in M_200
    X = p.X
    T_e = T_vir_calc(p, Mh * M_sun, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    nu = log(ustrip(X_to_nu(X)))
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    logMh = log10(Mh)
    dI = p.szpack_interp(t, nu)*u"MJy/sr"
    rsz_factor_I_over_y = (dI/(p.τ * θ_e))
    
    x₀ = cos(δ₀) * cos(α₀)
    y₀ = cos(δ₀) * sin(α₀) 
    z₀ = sin(δ₀)

    @inbounds for j in j_start:j_stop
        for i in i_start:i_stop
            x₁ = psa.cos_δ[i,j] * psa.cos_α[i,j]
            y₁ = psa.cos_δ[i,j] * psa.sin_α[i,j]
            z₁ = psa.sin_δ[i,j]
            d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
            θ = acos(clamp(1 - d² / 2, -one(T), one(T)))
            θ = max(θmin, θ)  # clamp to minimum θ
            y = exp(p.y_interp(log(θ), z, logMh))
            m[i,j] += (θ < θmax) * ustrip(u"MJy/sr", rsz_factor_I_over_y) * y
        end
    end
end


function profile_paint_szp!(m::HealpixMap{T, RingOrder}, p::Battaglia16SZPackProfile, 
            α₀, δ₀, w::HealpixProfileWorkspace, z, Mh, θmax) where T
    ϕ₀ = α₀
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    θmin = max(exp(first(first(interp_model.itp.ranges))), w.θmin)

    X = p.X
    T_e = T_vir_calc(p, Mh * M_sun, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    nu = log(ustrip(X_to_nu(X)))
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    logMh = log10(Mh)
    dI = p.szpack_interp(t, nu)*u"MJy/sr"
    rsz_factor_I_over_y = (dI/(p.τ * θ_e))
        
    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(clamp(1 - d² / 2, -one(T), one(T)))
        θ = max(θmin, θ)  # clamp to minimum θ
        y = exp(p.y_interp(log(θ), z, logMh))
        m.pixels[ir] += (θ < θmax) * ustrip(u"MJy/sr", rsz_factor_I_over_y) * y
    end
end


function paint_szp!(m::HealpixMap{T, RingOrder}, p::Battaglia16SZPackProfile, ws::Vector{W}, masses::AV, 
                    redshifts::AV, αs::AV, δs::AV) where {T, W <: HealpixProfileWorkspace, AV}
    m .= 0.0

    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_szp!(m, p, ws[i], masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_szp!(m, p, ws[i], masses, redshifts, αs, δs, i1:i2)
    end
end



# multi-halo painting utilities
function paint_szp!(m, p::XGPaint.AbstractProfile, psa, 
                masses::AV, redshifts::AV, αs::AV, δs::AV, irange::AbstractUnitRange) where AV
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax_ = θmax(p, mh * XGPaint.M_sun, z)
        profile_paint_szp!(m, p, α₀, δ₀, psa, z, mh, θmax_)
        
    end
end

function paint!(m, p::Battaglia16SZPackProfile, workspace, masses::AV, 
                        redshifts::AV, αs::AV, δs::AV)  where AV
    fill!(m, 0)
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    if N_sources < 2Threads.nthreads()  # don't thread if there are not many sources
        return paint_szp!(m, p, workspace, masses, redshifts, αs, δs, 1:N_sources)
    end
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_szp!(m, p, workspace, masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_szp!(m, p, workspace, masses, redshifts, αs, δs, i1:i2)
    end
end
