module XGPaint

using Distributed
using HDF5
using SharedArrays
using Healpix
using PoissonRandom
using Interpolations
using QuadGK
using Base.GC
using Roots
using Cosmology
using Unitful
using UnitfulAstro
using Random

# set different seeds for worker IDs
Random.seed!(myid() + trunc(Int64, time()))


"""
    CIBModel{T}(model parameters...)

Define CIB model parameters. Defaults are from Viero et al. 2013.

```@example
model = CIBModel{Float32}(shang_Mpeak=10^12.4)
```
"""
Base.@kwdef struct CIBModel{T}
    Inu_norm::T     = 0.3180384
    nside::T        = 4096
    min_redshift::T = 0.0
    max_redshift::T = 4.5
    min_mass::T     = 1e12
    box_size::T     = 40000
    hod::String     = "shang"
    LM::String      = "Planck2013"
    shang_zplat::T  = 2.0
    shang_Td::T     = 20.7
    shang_beta::T   = 1.6
    shang_eta::T    = 2.4
    shang_alpha::T  = 0.2
    shang_Mpeak::T  = 10^12.3
    shang_sigmaM::T = 0.3
end

function test(x::CIBModel{T}) where T
   return x.min_mass * 5
end


export CIBModel, test




end # module
