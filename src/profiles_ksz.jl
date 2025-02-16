

struct RKSZpackProfile{T,C,I1,I2,I3,I4} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
    X::T
    y_interp::I1
    tau_interp::I2
    szpack_interp_ksz::I3
    szpack_interp_T0::I4
end


function RKSZpackProfile(y_interp, tau_interp, szpack_interp_ksz, szpack_interp_T0;
        Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, x::T=2.6408) where T
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, Neff=3.046, OmegaM=OmegaM)
    @assert isangletypeparameter(tau_interp.model)
    return RKSZpackProfile(f_b, cosmo, x, y_interp, tau_interp, 
        szpack_interp_ksz, szpack_interp_T0)
end



"""
Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
"""
function SZpack_ksz(𝕡, M_200, z, r, vel; τ=0.01, mu = 1.0, showT=true)

    X = 𝕡.X
    T_e = T_vir_calc(𝕡, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)

    # use velocity magnitude to determine direction along line-of-sight
    if sign(vel) < 0
        mu *= -1
    end
    
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(X_to_nu(X)))
    vel = abs(uconvert(NoUnits, vel/constants.c_0))
    # need to take absolute value of v to make sure v is within bounds of interpolator
    
    # Term 1
    dI_1 = uconvert(u"kg*s^-2",(𝕡.szpack_interp_ksz(t, vel, mu, nu) * u"MJy/sr" - 
        𝕡.szpack_interp_T0(vel, mu, nu) * u"MJy/sr") / τ)
    y = compton_y(𝕡.y_interp.model, M_200, z, r)
    I_1 = uconvert(u"kg*s^-2",y * (dI_1/(θ_e)))
    
    # Term 2
    dI_2 = uconvert(u"kg*s^-2", (𝕡.szpack_interp_T0(vel, mu, nu) * u"MJy/sr")/τ)
    tau = XGPaint.tau(𝕡.tau_interp.model, r, M_200, z)
    I_2 = uconvert(u"kg*s^-2", dI_2 * tau)
    
    I = I_1 + I_2
    T = I/uconvert(u"kg*s^-2", abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X) / 
        (constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))

    if showT==true
        return abs(T)
    else
        return I
    end
end


function non_rel_ksz(𝕡, M_200, z, r, vel; mu = 1.0)
    """
    Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
    """

    # use velocity magnitude to determine direction along line-of-sight
    if vel < 0
        mu *= -1
    end
    
    vel = abs(ustrip(vel/uconvert(u"km/s",constants.c_0)))
    tau = XGPaint.tau(𝕡, r, M_200, z) #0.01 #XGPaint.tau_ksz(𝕡, M_200, z, r)

    # NON REL kSZ = tau * v/c (i.e. vel)
    T = tau*vel*mu

    return abs(T)
end


function profile_paint!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, p::RKSZpackProfile,
                        α₀, δ₀, psa::CarClenshawCurtisProfileWorkspace, 
                        z, Ms, vel, θmax) where T
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
            θ = max(θmin, θ)  # clamp to minimum θ
            y = exp(p.y_interp(log(θ), z, log10(Ms)))

            # m[i,j] += ifelse(θ < θmax, 
            #                  ,
            #                  zero(T))
        end
    end
end
