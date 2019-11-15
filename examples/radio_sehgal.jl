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

##
using DelimitedFiles
using Unitful
using UnitfulAstro
using Cosmology
using PyPlot
using StatsBase

##

figure(figsize=(8,8))

sehgal_red = readdlm("examples/data/sehgal_figure8_red.txt", ',', Float64, '\n')
sehgal_blue = readdlm("examples/data/sehgal_figure8_blue.txt", ',', Float64, '\n')
sehgal_green = readdlm("examples/data/sehgal_figure8_green.txt", ',', Float64, '\n')

vol = ustrip(u"Mpc^3",comoving_volume(cosmo, 0.3))
bins = (range(23, stop=29, step=0.5));
mids = (bins[2:end] .+ bins[1:end-1]) ./ 2.0
h = fit(Histogram, log10.(sources.L_I_151), bins)
plot(mids, h.weights ./ diff((bins)) ./ vol, label="websky draws I" )
h = fit(Histogram, log10.(sources.L_II_151), bins)
plot(mids, h.weights ./ diff((bins)) ./ vol, label="websky draws II" )

plot( log10.(sehgal_red[:,1]), sehgal_red[:,2] , "rp",
    label="fig 8. model \$\\times (P/P_0)\$", alpha=0.6)
plot( log10.(sehgal_green[:,1]), sehgal_green[:,2] , color="green", "^",
    label="fig. 8 FR I \$\\times (P/P_0)\$", alpha=0.6)
plot( log10.(sehgal_blue[:,1]), sehgal_blue[:,2] , color="blue", "s",
    label="fig. 8 FR II \$\\times (P/P_0)\$", alpha=0.6)

ylabel("counts / log bin / Mpc\$^3\$")
xlabel("log L [W/Hz]")
yscale("log")

legend()
gcf()

##
m = Map{Float64, RingOrder}(radio_model.nside)
result_map = m.pixels
nu_obs = 145f9
##

f = paint!(radio_model, sources, result_map, nu_obs)


##

##
# process centrals for this frequency
Threads.@threads for i in 1:sources.
    nu = (one(T) + sources.redshift_cen[i]) * nu_obs
    result_map[sources.hp_ind_cen[i]] += l2f(
        sources.lum_cen[i] * nu2theta(
            nu, sources.redshift_cen[i], model),
        sources.dist_cen[i], sources.redshift_cen[i])
end


    # set up indices for sources
    source_ind_I = generate_subhalo_indices(nsources_I)
    source_ind_II = generate_subhalo_indices(nsources_II)

##
"""
Paint a source catalog onto a map.
"""
function paint!(; nu_obs::T, result_map, sources, model::AbstractCIBModel) where T

    result_map .= zero(T)  # prepare the frequency map

    # process centrals for this frequency
    Threads.@threads for i in 1:sources.N_cen
        nu = (one(T) + sources.redshift_cen[i]) * nu_obs
        result_map[sources.hp_ind_cen[i]] += l2f(
            sources.lum_cen[i] * nu2theta(
                nu, sources.redshift_cen[i], model),
            sources.dist_cen[i], sources.redshift_cen[i])
    end

    # process satellites for this frequency
    Threads.@threads for i in 1:sources.N_sat
        nu = (one(T) + sources.redshift_sat[i]) * nu_obs
        result_map[sources.hp_ind_sat[i]] += l2f(
            sources.lum_sat[i] * nu2theta(
                nu, sources.redshift_sat[i], model),
            sources.dist_sat[i], sources.redshift_sat[i])
    end
end


##








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
