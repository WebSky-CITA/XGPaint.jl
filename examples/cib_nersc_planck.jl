using XGPaint
using Healpix
using Random
Random.seed!(3)

## Load halos from HDF5 files, establish a CIB model and cosmology
halo_pos, halo_mass = read_halo_catalog_hdf5(
    ENV["SCRATCH"] * "/websky_halos-light.hdf5")
# @load "/home/zequnl/testsamp.jld2" halo_mass halo_pos

cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = CIB_Planck2013{Float32}()

## Allocate some arrays and file them up for centrals and satellites
@time sources = generate_sources(model, cosmo, halo_pos, halo_mass);

## Deposit the sources into maps
fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)
m = HealpixMap{Float64,RingOrder}(model.nside)

for freq in [
    "18.7", "21.6",
    "24.5", "27.3", "30.0", "35.9",
    "41.7", "44.0", "47.4", "63.9", "67.8",
    "70.0", "73.7", "79.6", "90.2", "100", "111",
    "129", "143", "153", "164", "189", "210", "217",
    "232", "256", "275", "294", "306", "314", "340",
    "353", "375", "409", "467", "525", "545", "584",
    "643", "729", "817", "857", "906", "994", "1080"]

    @time begin
        XGPaint.paint!(m, parse(Float32, freq) * 1.0f9, model, sources,
            fluxes_cen, fluxes_sat)
        Healpix.saveToFITS(m, "!/global/cscratch1/sd/xzackli/cib/cib$(freq).fits")
    end
end
