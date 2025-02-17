
```@meta
CurrentModule = XGPaint
```

# The Sunyaev–Zeldovich Effect

Let's make some maps of the thermal Sunyaev-Zeldovich (SZ) effect. First, let's load up the example halo catalog included in this package. These will be automatically downloaded the first time you load them.

```@example tsz
using XGPaint, Plots

# example , ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
print("Number of halos: ", length(halo_mass))
```

This small catalog is limited to a relatively small patch of the sky. Before we generate some SZ maps, let's take a look at the halo mass distribution.

```@example tsz
b = 10.0 .^ (11:0.25:16)
histogram(halo_mass, bins=b, xaxis=(:log10, (1e11, 1e16)), 
    yscale=:log10, label="", xlabel="Halo mass (solar masses)", ylabel="counts")
```

Also notice the steep dropoff in counts at ``10^{12} M_{\odot}``, corresponding to the halo resolution of the Websky simulation from which these halos were taken. Fortunately, the SZ effect is dominated by the most massive halos.

To allow for safe threaded painting on a single map, we'll sort the halo catalog by declinations.

```@example tsz
ra, dec, redshift, halo_mass = sort_halo_catalog(ra, dec, redshift, halo_mass);
```

This package relies on a precomputed interpolation table for performance when making large sky maps of the SZ effect. For this tutorial, we'll load a profile and the associated interpolator from disk.

```@example tsz
y_model_interp = XGPaint.load_precomputed_battaglia()
print(y_model_interp)
```

!!! note

    If you want to generate your own model (such as varying cosmology), you would instead
    configure the appropriate model and then build an interpolator. The interpolator generation is multithreaded and takes about five minutes on 16 cores. 

    ```
    y_model = Battaglia16ThermalSZProfile(Omega_c=0.2589, Omega_b=0.0486, h=0.6774)
    y_model_interp = build_interpolator(y_model, cache_file="cached_b16.jld2", overwrite=true)
    ```
    This will save the results to the cache_file as well. If you want to load a result from disk, you can specify `overwrite=false`[^1].

[^1]: 
     For maps with many pixels and halos, building the interpolator is only about 10% of the map generation cost. If you have a use case where you instead need to vary cosmology on very small patches, contact the developers!


Now, we'll set up the map to generate. We will construct a small CAR (Clenshaw-Curtis variant) patch on the sky. You have to construct a workspace **for each new sky patch**.

```@example tsz
using Pixell
box = [4.5   -4.5;           # RA
       -3     3] * Pixell.degree  # DEC
shape, wcs = geometry(CarClenshawCurtis{Float64}, box, 0.5 * Pixell.arcminute)
m = Enmap(zeros(shape), wcs)

# construct on workspace on this patch
workspace = profileworkspace(shape, wcs)
```

Now it's time to apply the model to the map.
```@example tsz
@time paint!(m, workspace, y_model_interp, halo_mass, redshift, ra, dec)
plot(log10.(m), c = :thermal)
```


## Healpix 

To generate Healpix maps, you'll need the [Healpix.jl](https://github.com/ziotom78/Healpix.jl) package.

```@example tsz
using Healpix

nside = 4096
m_hp = HealpixMap{Float64,RingOrder}(nside)
max_radius = deg2rad(5.0)  # maximum radius to consider for the profile
w = HealpixProfileWorkspace(nside, max_radius)
@time paint!(m_hp, w, y_model_interp, halo_mass, redshift, ra, dec)
# Healpix.saveToFITS(m_hp, "!y.fits", typechar="D")  # to save, uncomment

plot(m_hp)
```
