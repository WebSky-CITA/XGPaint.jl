# XGPaint.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xzackli.github.io/XGPaint.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xzackli.github.io/XGPaint.jl/dev)
[![Build Status](https://github.com/xzackli/XGPaint.jl/workflows/CI/badge.svg)](https://github.com/xzackli/XGPaint.jl/actions)
[![codecov](https://codecov.io/gh/xzackli/XGPaint.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xzackli/XGPaint.jl)
<!-- [![Coveralls](https://coveralls.io/repos/github/xzackli/XGPaint.jl/badge.svg?branch=master)](https://coveralls.io/github/xzackli/XGPaint.jl?branch=master) -->


`XGPaint.jl` paints maps of extragalactic foregrounds using halo catalogs.

## Example (Planck 2013 CIB)

```julia
using XGPaint
using Healpix

## Load halos from HDF5 files, establish a CIB model and cosmology
halo_pos, halo_mass = read_halo_catalog_hdf5(
    "/home/zequnl/websky_halos-light.hdf5")
cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = CIBModel_Planck2013{Float32}()

## Allocate some arrays and file them up for centrals and satellites
sources = generate_sources(model, cosmo, halo_pos, halo_mass);

## Deposit the sources into maps
m = Map{Float64, RingOrder}(model.nside)
@time result = XGPaint.paint!(nu_obs=143.0f9,
    result_map=m.pixels, sources=sources, model=model)
Healpix.saveToFITS(m, "/media/data/cib143.fits")
```
