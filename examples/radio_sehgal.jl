using XGPaint
using Healpix
using JLD2


## Load halos from HDF5 files, establish a CIB model and cosmology
cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
radio_model = Radio_Sehgal2009{Float32}()

# halo_pos, halo_mass = read_halo_catalog_hdf5(
#     "/home/zequnl/websky_halos-light.hdf5")

@load "examples/data/lowz.jld2" halo_pos halo_mass


##
sources = generate_sources(radio_model, cosmo, halo_pos, halo_mass);

##
sources.L_I_151
##


using DelimitedFiles
using Unitful
using UnitfulAstro
using Cosmology
using PyPlot
using StatsBase

##

sehgal_red = readdlm("examples/data/sehgal_figure8_red.txt", ',', Float64, '\n')
sehgal_blue = readdlm("examples/data/sehgal_figure8_blue.txt", ',', Float64, '\n')
sehgal_green = readdlm("examples/data/sehgal_figure8_green.txt", ',', Float64, '\n')

vol = ustrip(u"Mpc^3",comoving_volume(cosmo, 0.3))
bins = (range(23, stop=29, step=0.5));
mids = (bins[2:end] .+ bins[1:end-1]) ./ 2.0


figure(figsize=(8,8))

h = fit(Histogram, log10.(sources.L_I_151), bins)
plot(mids, h.weights ./ diff((bins)) ./ vol, label="websky draws I" )
h = fit(Histogram, log10.(sources.L_II_151), bins)
plot(mids, h.weights ./ diff((bins)) ./ vol, label="websky draws II" )

plot( log10.(sehgal_red[:,1]), sehgal_red[:,2] , "rp", label="fig 8. model \$\\times (P/P_0)\$", alpha=0.6)
plot( log10.(sehgal_green[:,1]), sehgal_green[:,2] , color="green", "^", label="fig. 8 FR I \$\\times (P/P_0)\$", alpha=0.6)
plot( log10.(sehgal_blue[:,1]), sehgal_blue[:,2] , color="blue", "s", label="fig. 8 FR II \$\\times (P/P_0)\$", alpha=0.6)

ylabel("counts / log bin / Mpc\$^3\$")
xlabel("log L [W/Hz]")
yscale("log")

# plot( bins[2:end], 1e23 .* (10 .^bins[2:end]) .^ -1 )
legend()
gcf()
# ylim(2e-10, 1e-7)
# xticks(10.0 .^ collect(24:29))


# model = CIBModel_Planck2013{Float32}()
#
# ## Allocate some arrays and file them up for centrals and satellites
# sources = generate_sources(model, cosmo, halo_pos, halo_mass);
#
# ## Deposit the sources into maps
# m = Map{Float64, RingOrder}(model.nside)
# @time result = XGPaint.paint!(nu_obs=143.0f9,
#     result_map=m.pixels, sources=sources, model=model)
# Healpix.saveToFITS(m, "/media/data/cib143.fits")
