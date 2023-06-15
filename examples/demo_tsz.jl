using XGPaint, Plots

# example , ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
print("Number of halos: ", length(halo_mass))

##
b = 10.0 .^ (12:0.25:15)
histogram(halo_mass, bins=b, xaxis=(:log10, (1e12, 1e15)), yscale=:log10, label="", 
    xlabel="Halo mass (solar masses)", ylabel="number")

##
print("Precomputing the model profile grid.\n")

model = Battaglia16ThermalSZProfile(Omega_c=0.2589, Omega_b=0.0486, h=0.6774)
interp = build_interpolator(model, cache_file="cached_b16.jld2", overwrite=true)


ra, dec, redshift, halo_mass = sort_halo_catalog(ra, dec, redshift, halo_mass)

##
using Pixell
box = [4.5   -4.5;           # RA
       -3     3] * Pixell.degree  # DEC
shape, wcs = geometry(CarClenshawCurtis{Float64}, box, 0.5 * Pixell.arcminute)
m = Enmap(zeros(shape), wcs)

workspace = profileworkspace(shape, wcs)

@time paint!(m, model, workspace, interp, halo_mass, redshift, ra, dec)
# @time paint_map!(m, p, psa, sitp, halo_mass, redshift, ra, dec, 1:length(halo_mass))

##
plot(log10.(m), c = :thermal)

##
plot((m), clim=(-1e-5, 1e-5))



