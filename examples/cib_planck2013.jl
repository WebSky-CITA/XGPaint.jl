using XGPaint
using Healpix

## Load halos from HDF5 files, establish a CIB model and cosmology
halo_pos, halo_mass = read_halo_catalog_hdf5(
    "/tigress/zequnl/xgpaint/websky_halos-light.hdf5")
# @load "/home/zequnl/testsamp.jld2" halo_mass halo_pos

cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = CIB_Planck2013{Float32}()

## Allocate some arrays and file them up for centrals and satellites
@time sources = generate_sources(model, cosmo, halo_pos, halo_mass);

## Deposit the sources into maps
fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)
m = HealpixMap{Float64,RingOrder}(model.nside)

for freq in ["030", "090", "148", "219", "277", "350",
        "143", "217", "353", "545", "857"]
    @time begin
        XGPaint.paint!(m, parse(Float32, freq) * 1.0f9, model, sources,
            fluxes_cen, fluxes_sat)
        Healpix.saveToFITS(m, "!/tigress/zequnl/xgpaint/jl/cib/cib$(freq).fits")
    end
end
