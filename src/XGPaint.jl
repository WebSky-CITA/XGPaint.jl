
module XGPaint

using Interpolations
using QuadGK
using Roots
using LoopVectorization
using Cosmology
using Unitful, UnitfulAstro
using Parameters
using Random
using Distances
using Healpix, Pixell

import ThreadsX
import Distributions

include("./pixell.jl")
include("./model.jl")
include("./util.jl")
include("profile.jl")

include("./cib.jl")
include("./radio.jl")

export get_cosmology, read_halo_catalog_hdf5
export Radio_Sehgal2009, CIB_Planck2013
export paint!, generate_sources

end # module
