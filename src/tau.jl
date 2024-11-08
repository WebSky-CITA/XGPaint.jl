

# ANG is a Value type, Val(true) if we want to use angular units, and 
# Val(false) if we want to use physical units

struct BattagliaTauProfile{T,C,ANG} <: AbstractGNFW{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
end


function BattagliaTauProfile(; Omega_c::T=0.2589, Omega_b::T=0.0486, h::T=0.6774, angle=true) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(T, h=h, Neff=3.046, OmegaM=OmegaM)
    return BattagliaTauProfile{T, typeof(cosmo), angle}(f_b, cosmo)
end

function get_params(::BattagliaTauProfile{T}, M_200, z) where T
	z₁ = z + 1
	m = M_200 / (1e14M_sun)
    P₀ = 4.e3 * m^0.29 * z₁^(-0.66)
	α =  0.88 * m^(-0.03) * z₁^0.19
	β = -3.83 * m^0.04 * z₁^(-0.025)
	xc = 0.5
    γ = -0.2
    return (xc=T(xc), α=T(α), β=T(β), γ=T(γ), P₀=T(P₀))
end


function ρ_crit_comoving_h⁻²(p, z)
    return  (ρ_crit(p, z) ) / (1+z)^3 / p.cosmo.h^2
end

function r200c_comoving(p, m200c, z)
    rho_crit = ρ_crit(p, z) / (1+z)^3 
    return cbrt(m200c/ (4π/3 * rho_crit * 200)) 
end


function object_size(𝕡::BattagliaTauProfile{T,C,true}, physical_size, z) where {T,C}
    d_A = angular_diameter_dist(𝕡.cosmo, z)
    phys_siz_unitless = T(ustrip(uconvert(unit(d_A), physical_size)))
    d_A_unitless = T(ustrip(d_A))
    return atan(phys_siz_unitless, d_A_unitless)
end

function object_size(𝕡::BattagliaTauProfile{T,C,false}, physical_size, z) where {T,C}
    return physical_size
end


# returns a density, which we can check against Msun/Mpc² 
function rho_2d(p::BattagliaTauProfile, r, m200c, z)
    par = get_params(p, m200c, z)
    r200c = R_Δ(p, m200c, z, 200)
    X = r / object_size(p, r200c, z)  # either ang/ang or phys/phys
    rho_crit = ρ_crit_comoving_h⁻²(p, z)  # need to sort this out, it's all in comoving...
    result = par.P₀ * XGPaint._nfw_profile_los_quadrature(X, par.xc, par.α, par.β, par.γ)

    return result * rho_crit * (r200c * (1+z))
end

function ne2d(p::BattagliaTauProfile, r, m200c, z)
    me = constants.ElectronMass
    mH = constants.ProtonMass
    mHe = 4mH
    xH = 0.76
    nH_ne = 2xH / (xH + 1)
    nHe_ne = (1 - xH)/(2 * (1 + xH))
    factor = (me + nH_ne*mH + nHe_ne*mHe) / p.cosmo.h^2

    result = rho_2d(p, r, m200c, z)  # (Msun/h) / (Mpc/h)^2
    return result / factor
end

# r is either a physical or angular radius. unitful does not do little h, so physical radius 
# is always in Mpc, and angular radius is always in radians.
function tau(p, r, m200c, z)
    return constants.ThomsonCrossSection * ne2d(p, r, m200c, z) 
end

function profile_grid(𝕡::BattagliaTauProfile{T}, logθs, redshifts, logMs) where T

    N_logθ, N_z, N_logM = length(logθs), length(redshifts), length(logMs)
    A = zeros(T, (N_logθ, N_z, N_logM))

    Threads.@threads :static for im in 1:N_logM
        logM = logMs[im]
        M = 10^(logM) * M_sun
        for (iz, z) in enumerate(redshifts)
            for iθ in 1:N_logθ
                θ = exp(logθs[iθ])
                τ = tau(𝕡, θ, M, z)
                A[iθ, iz, im] = max(zero(T), y)
            end
        end
    end

    return logθs, redshifts, logMs, A
end