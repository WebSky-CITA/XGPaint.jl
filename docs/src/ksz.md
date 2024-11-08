
```@meta
CurrentModule = XGPaint
```

# The Kinetic Sunyaev–Zeldovich Effect

Let's make some maps of the kinetic Sunyaev-Zeldovich (SZ) effect. First, let's load up the example halo catalog included in this package. These will be automatically downloaded the first time you load them.

```@example ksz
using XGPaint, Plots

# example , ra and dec in radians, halo mass in M200c (Msun)
ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
print("Number of halos: ", length(halo_mass))
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
p = BattagliaTauProfile(Omega_c=0.267, Omega_b=0.0493,  h=0.6712)
```

We can now compute the integrated electron density, 

```@example ksz
zz = 2.5
Mnew = 93218298413772.23 * (M_sun / p.cosmo.h)  # Mcrit
tc = ne2d(p, 3u"Mpc" / p.cosmo.h, Mnew, zz) / (p.cosmo.h / 1u"Mpc")^2 + 0
@test abs(1 - tc / 1.4217533501414173e+69) < 1e-3
```

