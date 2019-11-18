using XGPaint
using Healpix

halo_pos, halo_mass = read_halo_catalog_hdf5(
    "/home/zequnl/websky_halos-light.hdf5")
## Load halos from HDF5 files, establish a CIB model and cosmology
cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
radio_model = Radio_Sehgal2009{Float32}(a_0=-1)

##

@time begin
    sources = generate_sources(radio_model, cosmo, halo_pos, halo_mass);
end;

##
@time begin
    m = Map{Float64,RingOrder}(radio_model.nside)
    paint!(m, 143f9, radio_model, sources)
    Healpix.saveToFITS(m, "/media/data/radio143.fits")
end
##
