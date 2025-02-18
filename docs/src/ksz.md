
```@meta
CurrentModule = XGPaint
```

# The Kinetic Sunyaev–Zeldovich Effect

Let's make some maps of the kinetic Sunyaev-Zeldovich (SZ) effect. First, let's load up the example halo catalog included in this package. These will be automatically downloaded the first time you load them.

```@example ksz
using XGPaint, Plots, Unitful, UnitfulAstro

# example , ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
print("Number of halos: ", length(halo_mass))

# the current test halos do not include velocities, so we randomly generate some for this
proj_v_over_c = randn(eltype(ra), length(halo_mass)) / 1000
```

This small catalog is limited to a relatively small patch of the sky. Before we generate some SZ maps, let's take a look at the halo mass distribution.

```@example ksz
b = 10.0 .^ (11:0.25:16)
histogram(halo_mass, bins=b, xaxis=(:log10, (1e11, 1e16)), 
    yscale=:log10, label="", xlabel="Halo mass (solar masses)", ylabel="counts")
```

To allow for safe threaded painting on a single map, we'll sort the halo catalog by declinations.

```@example ksz
ra, dec, redshift, halo_mass = sort_halo_catalog(ra, dec, redshift, halo_mass);
```


Now, we'll set up the map to generate. We will construct a small CAR (Clenshaw-Curtis variant) patch on the sky. You have to construct a workspace **for each new sky patch**.

```@example ksz
using Pixell
box = [4.5   -4.5;           # RA
       -3     3] * Pixell.degree  # DEC
shape, wcs = geometry(CarClenshawCurtis{Float64}, box, 0.5 * Pixell.arcminute)
m = Enmap(zeros(shape), wcs)
```

Now let's set up an electron profile.
```@example ksz
model = BattagliaTauProfile(Omega_c=0.267, Omega_b=0.0493,  h=0.6712)
```

We can now compute the integrated electron density, 

```@example ksz
workspace = profileworkspace(shape, wcs)

# this only needs to be done once
# to cache: interp = build_interpolator(model, cache_file="cached_btau.jld2", overwrite=false)
model_interp = build_interpolator(model)

m = Enmap(zeros(shape), wcs)
paint!(m, workspace, model_interp, halo_mass, redshift, ra, dec, proj_v_over_c)
plot(log10.(abs.(m)), c = :thermal)
```


## Healpix 

To generate Healpix maps, you'll need the [Healpix.jl](https://github.com/ziotom78/Healpix.jl) package.

```@example ksz
using Healpix

nside = 2048
m_hp = HealpixMap{Float64,RingOrder}(nside)
max_radius = deg2rad(5.0)  # maximum radius to consider for the profile
w = HealpixProfileWorkspace(nside, max_radius)

@time paint!(m_hp, w, model_interp, halo_mass, redshift, ra, dec, proj_v_over_c)
# Healpix.saveToFITS(m_hp, "!y.fits", typechar="D")  # to save, uncomment

plot(m_hp)
```
