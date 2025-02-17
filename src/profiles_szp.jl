

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
    y_model_interp::I1
    szpack_interp::I2
    τ::T
end

function SZPackRSZProfile(y_model_interp, x::T; τ=0.01, Omega_c=0.2589, 
        Omega_b=0.0486, h=0.6774, table_filename=rsz_szpack_table_filename()) where T
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    X = x
    szpack_interp = read_szpack_table(table_filename)
    return SZPackRSZProfile(f_b, cosmo, X, y_model_interp, szpack_interp, τ)
end

"""
    SZpack(model_szp, θ, M_200, z; τ=0.01, showT=true)

Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
"""
function SZpack(model_szp, θ, M_200, z; τ=0.01, showT=true)

    X = model_szp.X
    T_e = T_vir_calc(model_szp, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)

    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(X_to_nu(X)))
    dI = uconvert(u"kg*s^-2", model_szp.szpack_interp(t, nu)*u"MJy/sr")
    
    y = compton_y(model_szp.y_model_interp, θ, M_200, z)
    I = uconvert(u"kg*s^-2",y * (dI/(τ * θ_e)))
    T = I/uconvert(u"kg*s^-2",abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X) / 
        (constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))

    if showT==true
        return abs(T)
    else
        return I
    end
end

(model_szp::SZPackRSZProfile)(θ, M, z; τ=0.01, showT=true) = SZpack(
    model_szp, θ, M, z, τ=τ, showT=showT)


function I_to_T_mult_factor(X)
    return 1/(abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X) / 
        (constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))
end

function profile_paint!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, 
                        workspace::CarClenshawCurtisProfileWorkspace, 
                        model::SZPackRSZProfile, Mh, z, α₀, δ₀, θmax) where T
    X = model.X
    nu = log(ustrip(X_to_nu(X)))
    T_e = T_vir_calc(p, Mh * M_sun, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    dI = model.szpack_interp(t, nu) * u"MJy/sr"
    rsz_factor_I_over_y = dI / (model.τ * θ_e)
    profile_paint_generic!(m, workspace, model, Mh, z, α₀, δ₀, θmax, rsz_factor_I_over_y)
end

# like the usual paint, but use the sign of the null as the sign of the perturbation
function profile_paint!(m::HealpixMap{T, RingOrder}, 
        workspace::HealpixSerialProfileWorkspace, model::SZPackRSZProfile, 
        Mh, z, α₀, δ₀, θmax) where T
    X = model.X
    nu = log(ustrip(X_to_nu(X)))
    T_e = T_vir_calc(p, Mh * M_sun, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    nu = log(ustrip(X_to_nu(X)))
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    dI = model.szpack_interp(t, nu)*u"MJy/sr"
    rsz_factor_I_over_y = dI / (model.τ * θ_e)
    profile_paint_generic!(m, workspace, model, Mh, z, α₀, δ₀, θmax, rsz_factor_I_over_y)
end
