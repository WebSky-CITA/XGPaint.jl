using XGPaint

freqs = [143.0f9]
cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = CIBModel{Float32}()
##
halo_pos_serial, halo_mass_serial = read_halo_catalog_hdf5(
    "/home/zequnl/websky_halos-light.hdf5")
##
result = paint(freqs, model, cosmo, halo_pos_serial, halo_mass_serial)

##
using Healpix
m = Map{Float32, RingOrder}(model.nside)
m.pixels .= result;
Healpix.saveToFITS(m, "/media/data/cib$(freqs[1]).fits")
