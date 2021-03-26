
# Cosmic Infrared Background (CIB) 

We provide the Planck 2013 CIB model. The following code is a little more verbose than typical Julia code, as one has to repeatedly specify the type `Float32` when creating objects. This allows one to more easily fit the entire source catalog into memory.

## Sources

One first loads the halo positions and masses into memory. This package takes halo positions in the shape ``(3, N_{\mathrm{halos}})``, where the first dimension is the Cartesian coordinates ``x, y, z``.

```julia
using XGPaint
using Healpix

## Load halos from HDF5 files, establish a CIB model and cosmology
websky_directory = "/global/cfs/cdirs/sobs/www/users/Radio_WebSky"
halo_pos, halo_mass = read_halo_catalog_hdf5(
    "$(websky_directory)/websky_halos-light.hdf5")
```

Now one specifes the background cosmology and the source model [`CIB_Planck2013`](@ref).

```julia
# configuration objects
cosmo = get_cosmology(Float32; h=0.7, OmegaM=0.25)
model = CIB_Planck2013{Float32}()

# generate sources (healpix pixel, luminosities, etc. 
@time sources = generate_sources(model, cosmo, halo_pos, halo_mass);
```

This `sources` is a [NamedTuple](https://docs.julialang.org/en/v1/manual/types/#Named-Tuple-Types) with arrays for centrals,
* `hp_ind_cen`: healpix index of the central
* `lum_cen`: luminosity of the central
* `redshift_cen`: redshift of the central
* `dist_cen`: distance to the central

There are additionally arrays for the satellites,
* `hp_ind_sat`: healpix index of the satellite
* `lum_sat`: luminosity of the satellite
* `redshift_sat`: redshift of the satellite
* `dist_sat`: distance to the satellite

There are also two integers in `sources`  for the total number of centrals `N_cen` and total number of satellites `N_sat`.

## Map-making

Once these sources are generated, one needs to create some buffer arrays for map-making. The fluxes of the centrals and satellites are deposited into these arrays, before the map is generated.

```julia
# Create some empty arrays for the fluxes to be deposited
fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)
m = Map{Float64,RingOrder}(model.nside)  # create a Healpix map
```

These arrays are used by `paint!` to create maps. We then save those to disk. We add a time macro to get some info on how long it takes. Note that with NERSC's notoriously slow filesystem, writing to disk can take as long as generating the maps!

```julia
for freq in ["100", "143", "217" "353", "545"]
    @time paint!(m, parse(Float32, freq) * 1.0f9, model, 
        sources, fluxes_cen, fluxes_sat)
    saveToFITS(m, "!/global/cscratch1/sd/xzackli/cib/cib$(freq).fits")
end
```

## API
```@docs
CIB_Planck2013
generate_sources(
        model::XGPaint.AbstractCIBModel{T}, cosmo::Cosmology.FlatLCDM{T},
        halo_pos_inp::AbstractArray{TH,2}, halo_mass_inp::AbstractArray{TH,1};
        verbose=true) where {T, TH}
paint!(result_map::Healpix.Map{T_map, RingOrder},
        nu_obs, model::XGPaint.AbstractCIBModel{T}, sources,
        fluxes_cen::AbstractArray, fluxes_sat::AbstractArray) where {T_map, T}
```
