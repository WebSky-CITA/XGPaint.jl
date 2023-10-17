using XGPaint
using Healpix

# example , ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()

# sort
ra, dec, redshift, halo_mass = sort_halo_catalog(ra, dec, redshift, halo_mass)

print("Number of halos: ", length(halo_mass))

##
cosmo = get_cosmology(h=0.6774f0, OmegaM=0.3075f0)

x, y, z = XGPaint.ra_dec_redshift_to_xyz(ra, dec, redshift, cosmo)
halo_pos = [x'; y'; z';]

model = CIB_Planck2013{Float32}()

## Allocate some arrays and file them up for centrals and satellites
@time sources = generate_sources(model, cosmo, halo_pos, halo_mass);

## Deposit the sources into maps
fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)

using Pixell
box = [4.5   -4.5;           # RA
       -3     3] * Pixell.degree  # DEC
shape, wcs = geometry(CarClenshawCurtis{Float64}, box, 0.5 * Pixell.arcminute)
m = Enmap(zeros(Float32, shape), wcs)
XGPaint.paint!(m, 143.0f0 * 1.0f9, model, sources, fluxes_cen, fluxes_sat)

##
using Plots
plot(log10.(m), c=:coolwarm)




##


##
function get_basic_halo_properties(halo_pos::Array{T,2}, 
                                   cosmo) where T
    N_halos = size(halo_pos, 2)
    redshift = Array{T}(undef, N_halos)
    dist = Array{T}(undef, N_halos)

    r2z = XGPaint.build_r2z_interpolator(
        0.0, 4.5, cosmo)
    Threads.@threads :static for i in 1:N_halos
        dist[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
        redshift[i] = r2z(dist[i])
    end

    return dist, redshift
end

d̂, ẑ = get_basic_halo_properties([x'; y'; z'], cosmo)

θ, ϕ = XGPaint.get_angles([x'; y'; z'])
ra_ = ϕ
dec_ = π/2 .- θ

##
maximum((redshift .- ẑ)./redshift)


##
clip(x) = (x + 2π) % 2π

maximum(abs.((clip.(ra_) .- clip.(ra)))), maximum(abs.(dec_ .- dec))

##


##
