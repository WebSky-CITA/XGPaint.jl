

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

nu_to_X(nu) = (constants.h*nu)/(constants.k_B*T_cmb) + 0

struct SZPackRSZProfile{T,C,I1,I2} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
    X::T  # X = 2.6408 corresponding to frequency 150 GHz
    y_model_interp::I1
    szpack_interp::I2
    szpack_fiducial_tau::T
end

# forward min-theta to the y interpolator
compute_θmin(model::SZPackRSZProfile) = min(
    exp(first(first(model.y_model_interp.itp.ranges))))


function SZPackRSZProfile(y_model_interp, x::T; szpack_fiducial_tau=0.01, Omega_c=0.2589, 
        Omega_b=0.0486, h=0.6774, table_filename=rsz_szpack_table_filename()) where T
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, OmegaM=OmegaM)
    X = x
    szpack_interp = read_szpack_table(table_filename)
    return SZPackRSZProfile(f_b, cosmo, X, y_model_interp, szpack_interp, szpack_fiducial_tau)
end


"""
    SZpack(model_szp, θ, M_200, z; showT=false)

Outputs the integrated compton-y signal calculated using SZpack along the line of sight.
Note: M_200 requires units.
"""
function SZpack(model_szp, θ, M_200, z; showT=false)
    X = model_szp.X
    T_e = T_vir_calc(model_szp, M_200, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)

    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(uconvert(u"s^-1", X_to_nu(X))))
    dI = uconvert(u"kg*s^-2", model_szp.szpack_interp(t, nu)*u"MJy/sr")
    # y = compton_y(model_szp.y_model_interp.model, θ, M_200, z)
    Mh = M_200 / M_sun + 0
    y = model_szp.y_model_interp(θ, Mh, z)
    I = uconvert(u"kg*s^-2",y * (dI/(model_szp.szpack_fiducial_tau * θ_e)))
    T = I/uconvert(u"kg*s^-2",abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X) / 
        (constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))

    if showT==true
        return abs(T)
    else
        return I
    end
end

(model_szp::SZPackRSZProfile)(θ, M, z; showT=false) = ustrip(uconvert(
    u"MJy/sr", SZpack(model_szp, θ, M * M_sun, z, showT=showT)))


function I_to_T_mult_factor(X)
    return 1/(abs((2 * constants.h^2 * X_to_nu(X)^4 * ℯ^X) / 
        (constants.k_B * constants.c_0^2 * T_cmb * (ℯ^X - 1)^2)))
end

function compute_rsz_factor_I_over_y(model_szp::SZPackRSZProfile, Mh, z)
    X = model_szp.X
    T_e = T_vir_calc(model_szp, Mh * M_sun, z)
    θ_e = (constants.k_B*T_e)/(constants.m_e*constants.c_0^2)
    t = ustrip(uconvert(u"keV",T_e * constants.k_B))
    nu = log(ustrip(uconvert(u"s^-1", X_to_nu(X))))
    dI = uconvert(u"kg*s^-2", model_szp.szpack_interp(t, nu)*u"MJy/sr")
    rsz_factor_I_over_y = ustrip(uconvert(u"MJy/sr", 
        dI/(model_szp.szpack_fiducial_tau * θ_e)))
    return rsz_factor_I_over_y
end

# need to specify one for each type of map to prevent ambiguity
function profile_paint!(m::Enmap{T, 2, Matrix{T}, CarClenshawCurtis{T}}, 
                        workspace::CarClenshawCurtisProfileWorkspace, 
                        model_szp::SZPackRSZProfile, Mh, z, α₀, δ₀, θmax) where T
    rsz_factor_I_over_y = compute_rsz_factor_I_over_y(model_szp, Mh, z)
    profile_paint_generic!(m, workspace, model_szp.y_model_interp, Mh, z, α₀, δ₀, θmax, 
        rsz_factor_I_over_y)
end

function profile_paint!(m::HealpixMap{T, RingOrder}, 
        workspace::HealpixRingProfileWorkspace{T}, model::SZPackRSZProfile, 
        Mh, z, α₀, δ₀, θmax) where T
    rsz_factor_I_over_y = compute_rsz_factor_I_over_y(model, Mh, z)
    profile_paint_generic!(m, workspace, model.y_model_interp, Mh, z, α₀, δ₀, θmax, rsz_factor_I_over_y)
end
