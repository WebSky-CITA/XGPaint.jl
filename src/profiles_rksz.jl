

struct RKSZpackProfile{T,C,I1,I2,I3,I4} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
    X::T
    y_model_interp::I1
    tau_model_interp::I2
    szpack_interp_ksz::I3
    szpack_interp_T0::I4
    szpack_fiducial_tau::T
end


function RKSZpackProfile(y_model_interp, tau_interp, szpack_interp_ksz, szpack_interp_T0;
        Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, x::T=2.6408, szpack_fiducial_tau=0.01) where T
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, Neff=3.046, OmegaM=OmegaM)
    return RKSZpackProfile(f_b, cosmo, x, y_model_interp, tau_interp, 
        szpack_interp_ksz, szpack_interp_T0, szpack_fiducial_tau)
end


"""
    SZpack_rksz(model, r, M_200, z, vel; τ=0.01, mu = 1.0)

Computes the relativistic kSZ signal using the SZpack tables.
"""
function SZpack_rksz(model, r, M_200, z, vel, mu)

    X = model.X
    T_e = T_vir_calc(model, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)

    # use velocity magnitude to determine direction along line-of-sight
    mu *= sign(vel)
    
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(X_to_nu(X)))
    vel = abs(uconvert(NoUnits, vel/constants.c_0))
    # need to take absolute value of v to make sure v is within bounds of interpolator
    
    # Term 1
    dI_1 = (model.szpack_interp_ksz(t, vel, mu, nu) - 
        model.szpack_interp_T0(vel, mu, nu)) * u"MJy/sr" / model.szpack_fiducial_tau
    y = compton_y(model.y_model_interp.model, r, M_200, z)
    I_1 = y * dI_1 / θ_e
    
    # Term 2
    dI_2 = (model.szpack_interp_T0(vel, mu, nu) * u"MJy/sr") / model.szpack_fiducial_tau
    tau = compute_tau(model.tau_model_interp.model, r, M_200, z)
    I_2 = dI_2 * tau
    
    I = uconvert(u"MJy/sr", I_1 + I_2)
    return I
end

function T_over_Tcmb_from_I(model::RKSZpackProfile, I)
    X = model.X
    return I / abs((2 * constants.h^2 * X_to_nu(X)^4 * exp(X)) / 
            (constants.k_B * constants.c_0^2 * T_cmb * expm1(X)^2)) + 0
end

(model::RKSZpackProfile)(r, M, z, vel, mu) = SZpack_rksz(
    model, r, M * M_sun, z, vel, mu)


# rtkSZ takes in a last argument which is a tuple of the LOS velocity amplitude, and mu
function paintrange!(irange::AbstractUnitRange, m, workspace, model::RKSZpackProfile, 
                     masses, redshifts, αs, δs, velp_and_mu)
    nu = log(ustrip(uconvert(u"s^-1", X_to_nu(model.X))))
    for i in irange
        mass_with_units = masses[i] * XGPaint.M_sun
        θmax = compute_θmax(model, mass_with_units, redshifts[i])
        T_e = T_vir_calc(model, mass_with_units, redshifts[i])
        t = ustrip(uconvert(u"keV", T_e * constants.k_B))
        θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2) + 0
        vel, mu = velp_and_mu[i]  # amplitude of LOS velocity, and mu

        # in My/sr, no output units
        dI_1 = (model.szpack_interp_ksz(t, vel, mu, nu) - 
            model.szpack_interp_T0(vel, mu, nu)) / model.szpack_fiducial_tau
        profile_paint!(m, workspace, model.y_model_interp, 
                masses[i], redshifts[i], αs[i], δs[i], θmax, dI_1 / θ_e)
        
        # in My/sr, no output units
        dI_2 = model.szpack_interp_T0(vel, mu, nu) / model.szpack_fiducial_tau
        profile_paint!(m, workspace, model.tau_model_interp, 
                masses[i], redshifts[i], αs[i], δs[i], θmax, dI_2)

    end
end

