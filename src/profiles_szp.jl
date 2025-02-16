

function read_szpack_table(filename)
    table = readdlm(filename)
    nu_vector = LinRange(log(5.680062373019096*1e9),log(852.0093559528645*1e9),3000)
    temp_vector = LinRange(1.0e-3,75.0,100)
    szpack_interp = scale(Interpolations.interpolate(table, BSpline(Cubic(Line(OnGrid())))),
        temp_vector, nu_vector)
    return szpack_interp
end

function X_to_nu(X)
    return (X*constants.k_B*T_cmb)/constants.h
end

nu_to_X(nu) = (constants.h*nu)/(constants.k_B*T_cmb)

struct SZPackRSZProfile{T,C,I1,I2} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
    X::T  # X = 2.6408 corresponding to frequency 150 GHz
    model_y_interp::I1
    szpack_interp::I2
    τ::T
end

function SZPackRSZProfile(model_y_interp, x::T, τ=0.01; Omega_c=0.2589, 
        Omega_b=0.0486, h=0.6774, table_filename=rsz_szpack_table_filename()) where T
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    X = x
    szpack_interp = read_szpack_table(table_filename)
    return SZPackRSZProfile(f_b, cosmo, X, model_y_interp, szpack_interp, τ)
end

"""
    SZpack(model, M_200, z, r; τ=0.01, showT=true)

Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
"""
function SZpack(model_szp, θ, z, M_200; τ=0.01, showT=true)

    X = model_szp.X
    T_e = T_vir_calc(model_szp, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)

    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(X_to_nu(X)))
    dI = uconvert(u"kg*s^-2", model_szp.szpack_interp(t, nu)*u"MJy/sr")
    
    y = compton_y(model_szp.model_y_interp, θ, z, M_200)
    I = uconvert(u"kg*s^-2",y * (dI/(τ * θ_e)))
    T = I/uconvert(u"kg*s^-2",abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X) / 
        (constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))

    if showT==true
        return abs(T)
    else
        return I
    end
end

(model_szp::SZPackRSZProfile)(θ, z, M) = SZpack(model_szp, θ, z, M)


function I_to_T_mult_factor(X)
    return 1/(abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X) / 
        (constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))
end

function profile_paint!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, 
                        model::SZPackRSZProfile, 
                        workspace::CarClenshawCurtisProfileWorkspace, α₀, δ₀, 
                        z, Mh, θmax) where T
    # get indices of the region to work on
    i1, j1 = sky2pix(m, α₀ - θmax, δ₀ - θmax)
    i2, j2 = sky2pix(m, α₀ + θmax, δ₀ + θmax)
    i_start = floor(Int, max(min(i1, i2), 1))
    i_stop = ceil(Int, min(max(i1, i2), size(m, 1)))
    j_start = floor(Int, max(min(j1, j2), 1))
    j_stop = ceil(Int, min(max(j1, j2), size(m, 2)))
    θmin = compute_θmin(model)

    # needs mass in M_200
    X = model.X
    T_e = T_vir_calc(p, Mh * M_sun, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    nu = log(ustrip(X_to_nu(X)))
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    dI = model.szpack_interp(t, nu)*u"MJy/sr"
    rsz_factor_I_over_y = (dI/(model.τ * θ_e))
    
    x₀ = cos(δ₀) * cos(α₀)
    y₀ = cos(δ₀) * sin(α₀) 
    z₀ = sin(δ₀)

    @inbounds for j in j_start:j_stop
        for i in i_start:i_stop
            x₁ = workspace.cos_δ[i,j] * workspace.cos_α[i,j]
            y₁ = workspace.cos_δ[i,j] * workspace.sin_α[i,j]
            z₁ = workspace.sin_δ[i,j]
            d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
            θ = acos(clamp(1 - d² / 2, -one(T), one(T)))
            θ = max(θmin, θ)  # clamp to minimum θ
            y = model_szp.model_y_interp(θ, z, Mh)
            m[i,j] += (θ < θmax) * ustrip(u"MJy/sr", rsz_factor_I_over_y) * y
        end
    end
end


function profile_paint!(m::HealpixMap{T, RingOrder}, model_szp::SZPackRSZProfile, 
            w::HealpixProfileWorkspace, α₀, δ₀, z, Mh, θmax) where T
    ϕ₀ = α₀
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, m.resolution, θ₀, ϕ₀, θmax)
    θmin = max(compute_θmin(model), w.θmin)

    X = model.X
    T_e = T_vir_calc(p, Mh * M_sun, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    nu = log(ustrip(X_to_nu(X)))
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    dI = model.szpack_interp(t, nu)*u"MJy/sr"
    rsz_factor_I_over_y = (dI/(model.τ * θ_e))
        
    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap.pixels[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(clamp(1 - d² / 2, -one(T), one(T)))
        θ = max(θmin, θ)  # clamp to minimum θ
        y = model_szp.model_y_interp(θ, z, Mh)
        m.pixels[ir] += (θ < θmax) * ustrip(u"MJy/sr", rsz_factor_I_over_y) * y
    end
end
