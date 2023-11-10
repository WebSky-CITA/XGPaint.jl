

function read_szpack_table(filename)
    table = readdlm(filename)
    nu_vector = LinRange(log(35.6888844460172*1e9),log(5353.33266690298*1e9),3000)
    temp_vector = LinRange(1.0e-3,30.0,100)
    szpack_interp = scale(Interpolations.interpolate(table, BSpline(Cubic(Line(OnGrid())))), (temp_vector), (nu_vector))
    return szpack_interp
end

function SZpack(𝕡, M_200, z, r, τ=0.01)
    """
    Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
    """
    X = 𝕡.X
    T_e = T_vir_calc(𝕡, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    ω = (X*constants.k_B*T_cmb)/constants.ħ
    
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(uconvert(u"Hz",ω)))
    #if t > 30
    ##    t = 30
    #    println(t)
    #    println(M_200)
    #    println(z)
    #    println(X)
    #end
    dI = szpack_interp(t, nu)*u"MJy/sr"
   
    y = XGPaint.compton_y_rsz(𝕡, M_200, z, r)
    I = y * (dI/(τ * θ_e)) * (2π)^4
    T = I/abs((2 * constants.h^2 * ω^4 * ℯ^X)/(constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2))

    return abs(T)
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



function profile_paint_szp!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, p,
                        α₀, δ₀, psa::CarClenshawCurtisProfileWorkspace, 
                        sitp, z, Ms, θmax) where T
    # get indices of the region to work on
    i1, j1 = sky2pix(m, α₀ - θmax, δ₀ - θmax)
    i2, j2 = sky2pix(m, α₀ + θmax, δ₀ + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))
    
    X_0 = calc_null(p, Ms*M_sun, z)
    X = p.X
    if X > X_0
        sign = 1
    else
        sign = -1
    end 
    
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
                                 sign * exp(sitp(log(θ), z, log10(Ms))),
                                   zero(T))
        end
    end
end


"""Helper function to build a tSZ interpolator"""
function build_interpolator_szp(model::AbstractGNFW; cache_file::String="", 
                            N_logθ=512, pad=256, overwrite=true, verbose=true)

    if overwrite || (isfile(cache_file) == false)
        verbose && print("Building new interpolator from model.\n")
        rft = RadialFourierTransform(n=N_logθ, pad=pad)
        logθ_min, logθ_max = log(minimum(rft.r)), log(maximum(rft.r))
        prof_logθs, prof_redshift, prof_logMs, prof_y = profile_grid_szp(model; 
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
    sitp = scale(itp, prof_logθs, prof_redshift, prof_logMs)
    return sitp
end


function profile_paint_szp!(m::HealpixMap{T, RingOrder}, p,
            α₀, δ₀, w::HealpixProfileWorkspace, z, Mh, θmax) where T
    ϕ₀ = α₀
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    sitp = w.profile_real_interp
    
   X_0 = calc_null(p, Mh, z)
   X = p.X
   if X > X_0
       sign = 1
   else
       sign = -1
   end
    
    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(1 - d² / 2)
        θ = max(w.θmin, θ)  # clamp to minimum θ
        m.pixels[ir] += ifelse(θ < θmax, 
                                   sign * exp(sitp(log(θ), z, log10(Mh))),
                                    zero(T))
    end
end


# for rectangular pixelizations

# multi-halo painting utilities
function paint_szp!(m, p::XGPaint.AbstractProfile, psa, sitp, 
                masses::AV, redshifts::AV, αs::AV, δs::AV, irange::AbstractUnitRange) where AV
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax_ = θmax(p, mh * XGPaint.M_sun, z)
        profile_paint_szp!(m, p, α₀, δ₀, psa, sitp, z, mh, θmax_)
    end
end

function paint_szp!(m, p::XGPaint.AbstractProfile, psa, sitp, masses::AV, 
                        redshifts::AV, αs::AV, δs::AV)  where AV
    fill!(m, 0)
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_szp!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_szp!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2)
    end
end
