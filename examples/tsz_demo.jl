using XGPaint, Plots

# example , ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
print("Number of halos: ", length(halo_mass))

# ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
print(size(halo_mass))


##
b = 10.0 .^ (12:0.25:15)
histogram(halo_mass, bins=b, xaxis=(:log10, (1e12, 1e15)), yscale=:log10, label="", 
    xlabel="Halo mass (solar masses)", ylabel="number")

##

##
print("Precomputing the model profile grid.\n")


# p = XGPaint.BattagliaProfile(Omega_c=0.2589, Omega_b=0.0486, h=0.6774)
# beam stuff (not used in this particular script)
# N_logθ = 512
# rft = RadialFourierTransform(n=N_logθ, pad=256)

# model_file::String = "cached_battaglia.jld2"
# if isfile(model_file)
#     print("Found cached Battaglia profile model. Loading from disk.\n")
#     model = load(model_file)
#     prof_logθs, prof_redshift, prof_logMs, prof_y = model["prof_logθs"], 
#         model["prof_redshift"], model["prof_logMs"], model["prof_y"]
# else
#     print("Didn't find a cached profile model. Computing and saving.\n")
#     logθ_min, logθ_max = log(minimum(rft.r)), log(maximum(rft.r))
#     @time prof_logθs, prof_redshift, prof_logMs, prof_y = profile_grid(p; 
#         N_logθ=N_logθ, logθ_min=logθ_min, logθ_max=logθ_max)
#     save(model_file, Dict("prof_logθs"=>prof_logθs, 
#         "prof_redshift"=>prof_redshift, "prof_logMs"=>prof_logMs, "prof_y"=>prof_y))
# end


# set up a profile to paint
model, interp = XGPaint.load_precomputed_battaglia()



##

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


##

