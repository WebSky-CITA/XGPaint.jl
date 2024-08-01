```@meta
CurrentModule = XGPaint
```

# XGPaint

[XGPaint](https://github.com/xzackli/XGPaint.jl) paints maps of extragalactic foregrounds using halo models.
We provide CIB and radio models. The CIB model is described in the Websky paper, [Stein et al. 2020](https://arxiv.org/abs/2001.08787). The radio model is from Li et al. 2021, a paper in prep for the Websky suite, and is derived from [Sehgal et al. 2009](https://arxiv.org/abs/0908.0540).

The general workflow is:
1. Read halos
2. Specify a background cosmology and source model
3. Use the model to put sources in halos
4. Put the sources on maps

