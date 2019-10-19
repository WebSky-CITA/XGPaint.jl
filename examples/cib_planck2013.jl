using XGPaint
using Healpix

## Load halos from HDF5 files, establish a CIB model and cosmology
halo_pos, halo_mass = read_halo_catalog_hdf5(
    "/home/zequnl/websky_halos-light.hdf5")
cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = CIB_Planck2013{Float32}()

## Allocate some arrays and file them up for centrals and satellites
sources = generate_sources(model, cosmo, halo_pos, halo_mass);

## Deposit the sources into maps
m = Map{Float64, RingOrder}(model.nside)
@time result = XGPaint.paint!(nu_obs=143.0f9,
    result_map=m.pixels, sources=sources, model=model)
Healpix.saveToFITS(m, "/media/data/cib143.fits")
