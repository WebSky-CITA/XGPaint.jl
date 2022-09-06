
import PhysicalConstants.CODATA2018 as constants
using Unitful
const M_sun = 1.98847e30u"kg"
const P_e_factor = constants.σ_e / (constants.m_e * constants.c_0^2)
using Cosmology
using QuadGK
using XGPaint


abstract type AbstractProfile{T} end

struct BattagliaProfile{T,C} <: AbstractProfile{T}
    f_b::T  # Omega_b / Omega_c = 0.0486 / 0.2589
    cosmo::C
end

function BattagliaProfile(Omega_c::T=0.2589, Omega_b::T=0.0486) where {T <: Real}
    OmegaM=Omega_b+Omega_c
    f_b = Omega_b / OmegaM
    cosmo = get_cosmology(Float64, h=0.6774, OmegaM=OmegaM)
    return BattagliaProfile(f_b, cosmo)
end


function ρ_crit(𝕡, z)
    H_z = H(𝕡.cosmo, z)
    return uconvert(u"kg/m^3", 3H_z^2 / (8π * constants.G))
end

function R_Δ(𝕡, M_Δ, z, Δ=200)
    return (M_Δ / (4π/3 * Δ * ρ_crit(𝕡, z)))^(1/3)
end

function angular_size(𝕡::AbstractProfile, physical_size, z)
    d_A = angular_diameter_dist(𝕡.cosmo, z)

    # convert both to the same units and strip units for atan
    phys_siz_unitless = ustrip(uconvert(unit(d_A), physical_size))
    d_A_unitless = ustrip(d_A)
    return atan(phys_siz_unitless, d_A_unitless)
end

function generalized_nfw(x, xc, α, β, γ)
    x̄ = x / xc
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end

function _generalized_scaled_nfw(x̄, α, β, γ)
    return x̄^γ * (1 + x̄^α)^((β - γ) / α)
end


function get_params(::BattagliaProfile{T}, M_200, z) where T
	z₁ = z + 1
	m = M_200 / (1e14M_sun)
	P₀ = 18.1 * m^0.154 * z₁^-0.758
	xc = 0.497 * m^-0.00865 * z₁^0.731
	β = 4.35 * m^0.0393 * z₁^0.415
	α = 1
    γ = -0.3
    β = γ - α * β  # Sigurd's conversion from Battaglia to standard NFW
    return (xc=T(xc), α=T(α), β=T(β), γ=T(γ), P₀=T(P₀))
end

_tsz_y₁(x, _a) = (x*(_a+1))^(1/(_a+1))
_tsz_x₁(y, _a) = y^(_a+1)/(_a+1)
function _tsz_profile_los_quadrature(x, xc, α, β, γ; zmax=1e5, _a=8, rtol=eps())
    x² = x^2
    integral, err = quadgk(y -> y^_a * generalized_nfw(√(_tsz_x₁(y,_a)^2 + x²), xc, α, β, γ),
                      0.0, _tsz_y₁(zmax,_a), rtol=rtol)
    return 2integral
end

function dimensionless_P_profile_los(𝕡::BattagliaProfile{T}, M_200, z, r) where T
    par = get_params(𝕡, M_200, z)
    R_200 = R_Δ(𝕡, M_200, z, 200)
    x = r / angular_size(𝕡, R_200, z)
    return par.P₀ * _tsz_profile_los_quadrature(x, par.xc, par.α, par.β, par.γ)
end

"""Line-of-sight integrated electron pressure"""
P_e_los(𝕡, M_200, z, r) = 0.5176 * P_th_los(𝕡, M_200, z, r)

"""Line-of-sight integrated thermal pressure"""
P_th_los(𝕡, M_200, z, r) = constants.G * M_200 * 200 * ρ_crit(𝕡, z) * 
    𝕡.f_b / 2 * dimensionless_P_profile_los(𝕡, M_200, z, r)

function compton_y(𝕡, M_200, z, r)
    return P_e_los(𝕡, M_200, z, r) * P_e_factor
end


##
using BenchmarkTools, Test
p = BattagliaProfile()
@test abs(1 - compton_y(p, 3e15M_sun, 1.0, 0.0) / 0.000791435242101093) < 1e-5
@test abs(1 - compton_y(p, 5e11M_sun, 0.5, 1e-5) / 2.038204776991399e-08) < 1e-5
@test abs(1 - compton_y(p, 5e11M_sun, 1.0, 2e-5) / 1.6761760790073166e-08) < 1e-5
@test abs(1 - compton_y(p, 1e12M_sun, 2.0, 3e-5) /  4.646234867480562e-08) < 1e-5


##
rr = LinRange(0., deg2rad(4/60), 100)
plot(rad2deg.(rr) * 60,
    [compton_y(p, 1e12M_sun, 2.0, r) for r in rr],
    ylim=(0., 1e-10))

##
function gaussbeamA(fwhm::T, lmax::Int; pol=false) where T
    Bl = Array{T,1}(undef, lmax+1)
    fwhm²_to_σ² = 1 / (8log(T(2)))  # constant
    σ² = fwhm²_to_σ² * fwhm^2

    if pol
        for l = 0:lmax
            Bl[l+1] = exp(-(l * (l+1) - 4) * σ² / 2)
        end
    else
        for l = 0:lmax
            Bl[l+1] = exp(-l * (l+1) * σ² / 2)
        end
    end
    return Bl
end


##
using Pixell
rft = RadialFourierTransform(n=256, pad=128)
ells = rft.l
# lbeam = exp.(-ells .* (ells .+ 1) .* deg2rad(0.01/60)^2)

fwhm = deg2rad(1)
fwhm²_to_σ² = 1 / (8log(2))  # constant
σ² = fwhm²_to_σ² * fwhm^2

lbeam = @. exp(-ells * (ells+1) * σ² / 2)
plot(ells, lbeam, xlim=(0,1000))

# plot!(gaussbeamA(fwhm, 1000), ls=:dash)

##
# rprofs = [compton_y(p, 1e12M_sun, 2.0, r) for r in rft.r]
rprofs = [compton_y(p, 1e12M_sun, 2.0, r) for r in rft.r]
# rprofs = @. (rft.r[1]/4) ./ (rft.r)
##
using StaticArrays
z_min = 1e-3
z_max = 5.0
logM_min = 10  # 1e10 Msun
logM_max = 15

logθ_min = log(rft.r[begin+rft.pad])
logθ_max = log(rft.r[end-rft.pad])


N_z, N_logM, N_logθ = 64, 256, 512

##
Ms = LinRange(logM_min, logM_max, N_logM)
logθs = LinRange(logθ_min, logθ_max, N_logθ)
redshifts = LinRange(z_min, z_max, N_z)

LI = LinearIndices((1:N_logθ, 1:N_z, 1:N_logM))
# ys = zeros(N_logθ * N_z * N_logM)
# xs = zeros(SVector{3, Float64}, length(ys))
A = zeros((N_logθ, N_z, N_logM))

Threads.@threads for im in 1:N_logM
    logM = Ms[im]
    M = 10^(logM) * M_sun
    for (iz, z) in enumerate(redshifts)
        norm = compton_y(p, M, z, 0.0)
        for iθ in 1:N_logθ
            ii = LI[iθ, iz, im]
            θ = exp(logθs[iθ])
            A[iθ, iz, im] = compton_y(p, M, z, θ) / norm
            # ys[ii] = compton_y(p, M, z, θ)
            # xs[ii] = SA[logM, z, logθs[iθ]]
        end
    end
end


##
using Interpolations

itp = interpolate(A, BSpline(Cubic(Line(OnGrid()))))
sitp = scale(itp, logθs, redshifts, Ms)

##
@btime $sitp(-5.0, 1.0, 14.0)


##

ref(th, z, m) =  compton_y(p, 10^m * M_sun, z, exp(th)) / compton_y(p, 10^m * M_sun, z, 0.)

θs = exp.(logθs)
plot()
# plot(logθs, [abs(sitp(th, redshifts[10], Ms[14]) - ref(th, redshifts[10], Ms[14])) for th in logθs])
plot!(logθs, [abs(sitp(th, 0.5, 14.0) - ref(th, 0.5, 14.0)) for th in logθs])
plot!(logθs, [abs(sitp(th, redshifts[60], 14.1) - ref(th, redshifts[60], 14.1)) for th in logθs])

##

plot(logθs, [(sitp(th, 1e-3, 13.0)) for th in logθs])
plot!(logθs, [(ref(th,  0.5, 13.0)) for th in logθs])
plot!(logθs, [(sitp(th, 1.0, 13.0)) for th in logθs])
plot!(logθs, [(sitp(th, 2.0, 13.0)) for th in logθs])
# plot!(xscale=:log10)


##
plot(θs, A[:,1,64])
plot!(θs, A[:,2,64])
plot!(θs, A[:,32,64])
plot!(θs, A[:,64,64])
# plot!(θs, A[:,128,64])
# plot!(xscale=:log10, yscale=:log10)
plot!(xscale=:log10)
# plot!(xlim=(0, 0.1))


##

using FastChebInterp

# lb = minimum.((x[1], x[2], x[3]))
# ub = maximum.((x[1], x[2], x[3]))

c = chebregression(xs, ys, (10,10,10))


##

lprofs = real2harm(rft, rprofs)
lprofs .*= (lbeam)
rprofs2 = harm2real(rft, reverse(lprofs))


##
plot(rft.r[(begin+rft.pad):(end-rft.pad)], rprofs[(begin+rft.pad):(end-rft.pad)])
plot!(rft.r[(begin+rft.pad):(end-rft.pad)], rprofs2[(begin+rft.pad):(end-rft.pad)], ls=:dash, 
# xscale=:log10,
xlim=(0., 4e-2), 
ylim=(0.0, 2rprofs2[rft.pad])
)
vline!([fwhm/2])
# vline!([sol], ls=:dash)
# vline!([fwhm])


##


##
using DataInterpolations
yy = LinearInterpolation(rprofs2[rft.pad:(end-rft.pad)], rft.r[rft.pad:(end-rft.pad)])
using NonlinearSolve

ff(u,p) = yy(u) - rprofs2[rft.pad]/2
probB = NonlinearProblem(ff, 0.005)
sol = solve(probB, NewtonRaphson(), tol=1e-9).u

# lprofs .*= exp.(-rft.l .* (rft.l .+ 1) .* deg2rad(1/60)^2)

##

##
using DataInterpolations
r, rprofs = Pixell.unpad(rft, rft.r, rprofs)

##
plot(r, rprofs, xscale=:log)

##
plot(rft.r[64:end], rprofs[64:end], xscale=:log10)

##




##



##

using FastChebInterp


xr = rand(10000) * 10 # 10000 uniform random points in [0, 10]
c = chebregression(xr, f.(xr), 0, 10, 200) 


##
g(x) = sin(x[1] + cos(x[2]))
lb, ub = [1,3], [2, 4] # lower and upper bounds of the domain, respectively
x = chebpoints((10,20), lb, ub)
c = chebinterp(g.(x), lb, ub)

##
using StaticArrays
@btime $c(x) setup=(x=SA[1.5, 3.5])

##
@btime log(x) setup=(x=rand())



##


##
