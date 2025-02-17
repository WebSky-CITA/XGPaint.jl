

struct RKSZpackProfile{T,C,I1,I2,I3,I4} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
    X::T
    y_model_interp::I1
    tau_interp::I2
    szpack_interp_ksz::I3
    szpack_interp_T0::I4
end


function RKSZpackProfile(y_model_interp, tau_interp, szpack_interp_ksz, szpack_interp_T0;
        Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, x::T=2.6408) where T
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, Neff=3.046, OmegaM=OmegaM)
    @assert isangletypeparameter(tau_interp.model)
    return RKSZpackProfile(f_b, cosmo, x, y_model_interp, tau_interp, 
        szpack_interp_ksz, szpack_interp_T0)
end



"""
Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
"""
function SZpack_rksz(model, r, M_200, z, vel; τ=0.01, mu = 1.0, showT=true)

    X = model.X
    T_e = T_vir_calc(model, M_200, z)
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
    dI_1 = uconvert(u"kg*s^-2",(model.szpack_interp_ksz(t, vel, mu, nu) * u"MJy/sr" - 
        model.szpack_interp_T0(vel, mu, nu) * u"MJy/sr") / τ)
    y = compton_y(model.y_model_interp.model, r, M_200, z)
    I_1 = uconvert(u"kg*s^-2",y * (dI_1/(θ_e)))
    
    # Term 2
    dI_2 = uconvert(u"kg*s^-2", (model.szpack_interp_T0(vel, mu, nu) * u"MJy/sr")/τ)
    tau = tau(model.tau_interp.model, r, M_200, z)
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

(model::RKSZpackProfile)(r, M, z, vel; τ=0.01, mu = 1.0, showT=true) = SZpack_rksz(
    model, r, M, z, vel, τ=τ, mu=mu, showT=showT)

function non_rel_ksz(model, r, M_200, z, vel; mu = 1.0)

    # use velocity magnitude to determine direction along line-of-sight
    if vel < 0
        mu *= -1
    end
    
    vel = abs(ustrip(vel/uconvert(u"km/s",constants.c_0)))
    tau = tau(model, r, M_200, z) #0.01 #XGPaint.tau_ksz(model, M_200, z, r)

    # NON REL kSZ = tau * v/c (i.e. vel)
    T = tau*vel*mu

    return abs(T)
end

