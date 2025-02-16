
```@meta
CurrentModule = XGPaint
```

# The Relativistic Sunyaev–Zeldovich Effect

An example from Kuhn et al. (in prep). Begin with the tSZ example but on a smaller tile,

```@example rsz
using Unitful, UnitfulAstro, XGPaint, Plots, Pixell
import PhysicalConstants.CODATA2018 as constants

M_200 = 1e15 * XGPaint.M_sun
z = 0.4
θ = 60*4.84814e-6

box = [-30   30;           # RA
       -30   30] * Pixell.arcminute  # DEC
shape, wcs = geometry(Pixell.CarClenshawCurtis, box, 1 * Pixell.arcminute)
workspace = profileworkspace(shape, wcs)
model_tsz, y_interp = XGPaint.load_precomputed_battaglia();
nothing # hide
```


```@example rsz
tsz_snap = Enmap(zeros(shape), wcs)
mass_in_Msun = ustrip(M_200 / XGPaint.M_sun)
paint!(tsz_snap, model_tsz, workspace, y_interp, [mass_in_Msun], [z], [0.0], [0.0])
plot(log10.(tsz_snap))
```


```@example rsz
tsz_snap = Enmap(zeros(shape), wcs)
mass_in_Msun = ustrip(M_200 / XGPaint.M_sun)
paint!(tsz_snap, model_tsz, workspace, y_interp, [mass_in_Msun], [z], [0.0], [0.0])
plot(log10.(tsz_snap))
```



```@example rsz
i_ctr, j_ctr = round.(Int, sky2pix(tsz_snap, 0.0, 0.0))
ras = [pix2sky(tsz_snap, i, j_ctr)[1] for i in 1:size(tsz_snap,1)]
plot(ras, log10.(tsz_snap[:,j_ctr]), label="tSZ", xlabel="RA [rad]", ylabel="log10 y")
```



```@example rsz
# generate a snapshot of a cluster at a specific frequency
function cluster_snapshot(nu, shape, wcs, model_tsz, y_interp, workspace)
    X = nu_to_X(uconvert(u"s^-1",nu))
    m = Enmap(zeros(shape), wcs);

    # construct a new SZpack object for each frequency; this loads a table from disk.
    # keep these around if you're looping over multiple maps at fixed frequency
    model_szp = Battaglia16SZPackProfile(model_tsz, y_interp, X)

    paint!(m, model_szp, workspace, [mass_in_Msun], [z], [0.0], [0.0])
    return m
end

snap = cluster_snapshot(500u"GHz", shape, wcs, model_tsz, y_interp, workspace)
plot(log10.(snap))
```


```@example rsz
freqs = LinRange(30.0, 800.0, 100) .* u"GHz"
snaps = [cluster_snapshot(nu, shape, wcs, model_tsz, y_interp, workspace) for nu in freqs]
I_rsz = [snap[i_ctr, j_ctr] * u"MJy/sr" for snap in snaps]
I_tsz = similar(I_rsz)

y0 = tsz_snap[i_ctr, j_ctr]
for (i, nu) in enumerate(freqs)
    X = nu_to_X(nu)
    I_tsz[i]  = uconvert(u"MJy/sr",(X*(ℯ^X + 1)/(ℯ^X - 1) - 4) * y0 * 
        abs((2 * constants.h^2 * nu^4 * ℯ^X) / (
            constants.k_B * constants.c_0^2 * XGPaint.T_cmb * (ℯ^X - 1)^2)))
end

plot(freqs, I_tsz, label="tSZ")
plot!(freqs, I_rsz, label="rSZ")
```
