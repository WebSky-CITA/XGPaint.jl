
```@meta
CurrentModule = XGPaint
```

# The Sunyaev–Zeldovich Effect

Let's make some maps of the Sunyaev-Zeldovich (SZ) effect. First, let's load up the example halo catalog included in this package. These will be automatically downloaded the first time you call the function `load_example_halos`.

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

Also notice the steep dropoff in counts at ``10^{12} M_{\odot}``, corresponding to the halo resolution of the Websky simulation from which these halos were taken.

To allow for safe threaded painting on a single map, we'll sort the halo catalog by declinations.

```@example tsz
ra, dec, redshift, halo_mass = sort_halo_catalog(ra, dec, redshift, halo_mass);
```

This package relies on a precomputed interpolation table for performance when making large sky maps of the SZ effect. For this tutorial, we'll load a profile and the associated interpolator from disk.

```@example tsz
model, interp = XGPaint.load_precomputed_battaglia()
print(model)
```

Now, we'll set up the map to paint. We will construct a standard small CAR (Clenshaw-Curtis variant) patch on the sky. You have to construct a workspace **for each new sky patch**.

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
@time paint!(m, model, workspace, interp, halo_mass, redshift, ra, dec)
plot(log10.(m), c = :thermal)
```
