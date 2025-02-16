
struct Battaglia16RelativisticSZProfile{T,C} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
    X::T  # X = 2.6408 corresponding to frequency 150 GHz
end

function Battaglia16RelativisticSZProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, x::T=2.6408) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    X = x
    return Battaglia16RelativisticSZProfile(f_b, cosmo, X)
end

function dimensionless_P_profile_los_rsz(model::Battaglia16RelativisticSZProfile{T}, M_200, z, r) where T
    par = get_params(model, M_200, z)
    R_200 = R_Δ(model, M_200, z, 200)
    x = r / angular_size(model, R_200, z)
    return par.P₀ * _nfw_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
end

"""Line-of-sight integrated electron pressure"""
P_e_los_rsz(model, M_200, z, r) = 0.5176 * P_th_los_rsz(model, M_200, z, r)

"""Line-of-sight integrated thermal pressure"""
P_th_los_rsz(model, M_200, z, r) = constants.G * M_200 * 200 * ρ_crit(model, z) * 
    model.f_b / 2 * dimensionless_P_profile_los_rsz(model, M_200, z, r)

function compton_y_rsz(model, M_200, z, r)
    return P_e_los_rsz(model, M_200, z, r) * P_e_factor
end


function T_vir_calc(model,M,z::T) where T
   """
   Calculates the virial temperature for a given halo using Wang et al. 2007.
   """
    µ = 0.6  #µ is the mean molecular weight -> used the primordial abundance
    if z >= 1
        d_c = T(178)
    else
        d_c = T(356/(1 + z))
    end
    T_vir = 4.8e-3 * (M/M_sun)^(2/3) * (1 + z) * (model.cosmo.Ω_m/0.3)^(1/3) * (d_c/178)^(1/3) * (µ/0.59) * u"K"
    return T_vir
end

function T_mass_calc(model,M,z::T; scale_type="Ty", sim_type="combination") where T
    """
   Calculates the temperature for a given halo using https://arxiv.org/pdf/2207.05834.pdf.
   """
    E_z = model.cosmo.Ω_m*(1 + z)^3 + model.cosmo.Ω_Λ
    par_dict_scale = Dict([("Ty",[1.426,0.566,0.024]),("Tm",[1.207,0.589,0.003]),("Tsl",[1.196,0.641,-0.048])])
    par_dict_sim = Dict([("combination",[1.426,0.566,0.024]),("bahamas",[2.690,0.323,0.023]),("the300",[2.294,0.350,0.013]),("magneticum",[2.789,0.379,0.030]),("tng",[2.154,0.365,0.032])])
    
    if scale_type=="Ty"
        params = par_dict_sim[sim_type]
    else
        params = par_dict_scale[scale_type]
    end
    
    T_e = E_z^(2/3) * params[1] * (M/(10^14 *M_sun))^(params[2] + params[3] * log10(M/(10^14 * M_sun))) * u"keV"
    
    return T_e  
end

function rSZ(model, M_200, z, r; T_scale="virial", sim_type="combination", showT=true)
    """
    Calculates the integrated relativistic compton-y signal along the line of sight.
    """
    #X = (constants.ħ*ω)/(constants.k_B*T_cmb) # omega is standard frequency in Hz
    X = model.X
    if T_scale=="virial"
        T_e = T_vir_calc(model, M_200, z)
    elseif typeof(T_scale)==String
        T_e = uconvert(u"K",(T_mass_calc(model, M_200, z, scale_type=T_scale, sim_type=sim_type)/constants.k_B))
    else
        T_e = T_scale
    end
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    ω = (X*constants.k_B*T_cmb)/constants.h

    Xt = X*coth(X/2)
    St = X/(sinh(X/2))

    Y_0 = -4 + Xt
    Y_1 = -10 + (47/2)*Xt -(42/5)*Xt^2 + (7/10)*Xt^3 + St^2*((-21/5)+(7/5)*Xt)
    Y_2 = (-15/2) + (1023/8)*Xt - (868/5)*Xt^2 + (329/5)*Xt^3 - (44/5)*Xt^4 + (11/30)*Xt^5 +
        St^2*((-434/5) + (658/5)*Xt - (242/5)*Xt^2 + (143/30)*Xt^3) + St^4*((-44/4) + (187/60)*Xt)
    Y_3 = (15/2) + (2505/8)*Xt - (7098/5)*Xt^2 + (14253/10)*Xt^3 - (18594/35)*Xt^4 + (12059/140)*Xt^5 -
        (128/21)*Xt^6 + (16/105)*Xt^7 + St^2*((-7098/10) + (14253/5)*Xt - (102267/35)*Xt^2 +
        (156767/140)*Xt^3 - (1216/7)*Xt^4 + (64/7)*Xt^5) + St^4*((-18594/35) + (205003/280)*Xt -
        (1920/7)*Xt^2 + (1024/35)*Xt^3) + St^6*((-544/21) + (992/105)*Xt)
    Y_4 = (-135/32) + (30375/128)*Xt - (62391/10)*Xt^2 + (614727/40)*Xt^3 - (124389/10)*Xt^4 +
        (355703/80)*Xt^5 - (16568/21)*Xt^6 + (7516/105)*Xt^7 - (22/7)*Xt^8 + (11/210)*Xt^9 +
        St^2*((-62391/20) + (614727/20)*Xt - (1368279/20)*Xt^2 + (4624139/80)*Xt^3 - (157396/7)*Xt^4 +
        (30064/7)*Xt^5 - (2717/7)*Xt^6 + (2761/210)*Xt^7) + St^4*((-124389/10) + (6046951/160)*Xt -
        (248520/7)*Xt^2 + (481024/35)*Xt^3 - (15972/7)*Xt^4 + (18689/140)*Xt^5) + St^6*((-70414/21) +
        (465992/105)*Xt - (11792/7)*Xt^2 + (19778/105)*Xt^3) + St^8*((-628/7) + (7601/210)*Xt)

    prefac = ((X*ℯ^X)/(ℯ^X-1))*θ_e*(Y_0+θ_e*Y_1+θ_e^2*Y_2+θ_e^3*Y_3+θ_e^4*Y_4)
    y = compton_y_rsz(model, M_200, z, r)
    n = prefac * (constants.m_e*constants.c_0^2)/(T_e*constants.k_B) * y
    I = (X^3/(ℯ^X-1)) * (2*(constants.k_B*T_cmb)^3)/((constants.h*constants.c_0)^2) * n 
    T = I/abs((2 * constants.h^2 * ω^4 * ℯ^X)/(constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2))

    if showT==true
        return abs(T)
    else
        return I
    end
end

function calc_null(model, M_200, z)
    T_e = T_vir_calc(model, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    X_0 = 3.830*(1 + 1.1674*θ_e - 0.8533*(θ_e^2))
    
    return X_0
end


function profile_grid_rsz(model::AbstractGNFW{T}; N_z=256, N_logM=256, N_logθ=512, z_min=1e-3, z_max=5.0, 
              logM_min=11, logM_max=15.7, logθ_min=-16.5, logθ_max=2.5) where T

    logθs = LinRange(logθ_min, logθ_max, N_logθ)
    redshifts = LinRange(z_min, z_max, N_z)
    logMs = LinRange(logM_min, logM_max, N_logM)

    return profile_grid_rsz(model, logθs, redshifts, logMs)
end


function profile_grid_rsz(model::AbstractGNFW{T}, logθs, redshifts, logMs) where T

    N_logθ, N_z, N_logM = length(logθs), length(redshifts), length(logMs)
    A = zeros(T, (N_logθ, N_z, N_logM))

    Threads.@threads for im in 1:N_logM
        logM = logMs[im]
        M = 10^(logM) * M_sun
        for (iz, z) in enumerate(redshifts)
            for iθ in 1:N_logθ
                θ = exp(logθs[iθ])
                rsz = rSZ(model, M, z, θ)
                A[iθ, iz, im] = max(zero(T), rsz)
            end
        end
    end

    return logθs, redshifts, logMs, A
end



function profile_paint_rsz!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, model,
                        α₀, δ₀, psa::CarClenshawCurtisProfileWorkspace, 
                        z, Mh, θmax) where T
    # get indices of the region to work on
    i1, j1 = sky2pix(m, α₀ - θmax, δ₀ - θmax)
    i2, j2 = sky2pix(m, α₀ + θmax, δ₀ + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))
    θmin = compute_θmin(model)
    
    X_0 = calc_null(p, Mh*M_sun, z)
    X = model.X
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
            θ =  acos(clamp(1 - d² / 2, -one(T), one(T)))
            θ = max(θmin, θ)  # clamp to minimum θ
            m[i,j] += ifelse(θ < θmax, sign * model(θ, z, Mh), zero(T))
        end
    end
end


function profile_paint_rsz!(m::HealpixMap{T, RingOrder}, model,
            α₀, δ₀, w::HealpixProfileWorkspace, z, Mh, θmax) where T
    ϕ₀ = α₀
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    θmin = max(compute_θmin(model), w.θmin)

    X_0 = calc_null(p, Mh, z)
    X = model.X
    if X > X_0
        sign = 1
    else
        sign = -1
    end

    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ =  acos(clamp(1 - d² / 2, -one(T), one(T)))
        θ = max(θmin, θ)  # clamp to minimum θ
        m.pixels[ir] += ifelse(θ < θmax, sign * model(θ, z, Mh), zero(T))
    end
end


# for rectangular pixelizations

# multi-halo painting utilities
function paint_rsz!(m, p::XGPaint.AbstractProfile, psa, 
                masses::AV, redshifts::AV, αs::AV, δs::AV, irange::AbstractUnitRange) where AV
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax_ = compute_θmax(p, mh * XGPaint.M_sun, z)
        profile_paint_rsz!(m, p, α₀, δ₀, psa, z, mh, θmax_)
    end
end

function paint_rsz!(m, p::XGPaint.AbstractProfile, psa, masses::AV, 
                        redshifts::AV, αs::AV, δs::AV)  where AV
    fill!(m, 0)
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_rsz!(m, p, psa, masses, redshifts, αs, δs, i1:i2)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_rsz!(m, p, psa, masses, redshifts, αs, δs, i1:i2)
    end
end
