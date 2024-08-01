
# Cosmic Infrared Background (CIB) 

We provide the Planck 2013 CIB model. The following code is a little more verbose than typical Julia code, as one has to repeatedly specify the type `Float32` when creating objects. This allows one to more easily fit the entire source catalog into memory.


## Tutorial


```@example cib
using XGPaint, Plots, Pixell

# example , ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()

# sort
ra, dec, redshift, halo_mass = sort_halo_catalog(ra, dec, redshift, halo_mass)

print("Number of halos: ", length(halo_mass))
```

Next, we'll generate a cosmology. Note how we use Float32 throughout.

```@example cib

cosmo = get_cosmology(h=0.6774f0, OmegaM=0.3075f0)
x, y, z = XGPaint.ra_dec_redshift_to_xyz(ra, dec, redshift, cosmo)
halo_pos = [x'; y'; z';]

model = CIB_Planck2013{Float32}()
```


```@example cib
@time sources = generate_sources(model, cosmo, halo_pos, halo_mass);
fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)

using Pixell
box = [4.5   -4.5;           # RA
       -3     3] * Pixell.degree  # DEC
shape, wcs = geometry(CarClenshawCurtis{Float64}, box, 0.5 * Pixell.arcminute)
m = Enmap(zeros(Float32, shape), wcs)
XGPaint.paint!(m, 143.0f0 * 1.0f9, model, sources, fluxes_cen, fluxes_sat)
plot(log10.(m), c=:coolwarm)
```

## Sources from HDF5

To work with the Websky halo catalogs, you can't just use the example catalogs! One first loads the halo positions and masses into memory with [`read_halo_catalog_hdf5`](@ref). This package uses halo positions in the shape ``(3, N_{\mathrm{halos}})``, where the first dimension is the Cartesian coordinates ``x, y, z``.

```julia
using XGPaint
using Healpix

## Load halos from HDF5 files, establish a CIB model and cosmology
websky_directory = "/global/cfs/cdirs/sobs/www/users/Radio_WebSky"
halo_pos, halo_mass = read_halo_catalog_hdf5(
    "$(websky_directory)/websky_halos-light.hdf5")
```

Now one specifes the background cosmology with [`get_cosmology`](@ref) and the source model [`CIB_Planck2013`](@ref).

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

## HealpixMap-making

Once these sources are generated, one needs to create some buffer arrays for map-making. The fluxes of the centrals and satellites are deposited into these arrays, before the map is generated.

```julia
# Create some empty arrays for the fluxes to be deposited
fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)
m = HealpixMap{Float64,RingOrder}(model.nside)  # create a Healpix map
```

These arrays are used by `paint!` to create maps. We then save those to disk. We add a time macro to get some info on how long it takes. Note that with NERSC's notoriously slow filesystem, writing to disk can take as long as generating the maps!

```julia
for freq in ["100", "143", "217" "353", "545"]
    @time paint!(m, parse(Float32, freq) * 1.0f9, model, 
        sources, fluxes_cen, fluxes_sat)
    saveToFITS(m, "!/global/cscratch1/sd/xzackli/cib/cib$(freq).fits")
end
```

## Custom Models

You can make changes while reusing the XGPaint infrastructure by using Julia's multiple dispatch.
Create a custom type that inherits from `AbstractCIBModel`.
You must `import` the function you want to replace,
and then write your own version of the function which dispatches on your custom type.

```julia
 # import AbstractCIBModel and the functions you want to replace
import XGPaint: AbstractCIBModel, shang_z_evo 
using Parameters, Cosmology

# write your own type that is a subtype of AbstractCIBModel
@with_kw struct CustomCIB{T<:Real} <: AbstractCIBModel{T} @deftype T
    nside::Int64    = 4096
    hod::String     = "shang"
    Inu_norm     = 0.3180384
    min_redshift = 0.0
    max_redshift = 5.0
    min_mass     = 1e12
    box_size     = 40000

    # shang HOD
    shang_zplat  = 2.0
    shang_Td     = 20.7
    shang_beta   = 1.6
    shang_eta    = 2.4
    shang_alpha  = 0.2
    shang_Mpeak  = 10^12.3
    shang_sigmaM = 0.3
    shang_Msmin  = 1e11
    shang_Mmin   = 1e10
    shang_I0     = 46

    # jiang
    jiang_gamma_1    = 0.13
    jiang_alpha_1    = -0.83
    jiang_gamma_2    = 1.33
    jiang_alpha_2    = -0.02
    jiang_beta_2     = 5.67
    jiang_zeta       = 1.19
end

# dispatch on custom type CustomCIB. this particular change sets L=0
function shang_z_evo(z::T, model::CustomCIB) where T
    return zero(T)
end

# make a cosmology and our custom source model
cosmo = get_cosmology(Float32; h=0.7, OmegaM=0.3)
custom_model = CustomCIB{Float32}()

# for this test, do it only on a subset of the halos
sources = generate_sources(custom_model, cosmo, halo_pos[:,1:10], halo_mass[1:10])

print(sources.lum_cen)
```

This particular change zeros out the luminosities, and indeed you should see the result is an array of zeroes.
