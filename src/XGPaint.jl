
module XGPaint

using Interpolations
using QuadGK
using Roots
using Cosmology
using Unitful, UnitfulAstro
using Parameters
using Random
using Healpix
using PhysicalConstants
using Pixell
using Healpix: checkPixelRing
using SpecialFunctions
using JLD2, FileIO

import ThreadsX
import Distributions
import DataInterpolations  # used for uneven spaced interpolators

include("./util.jl")
include("./model.jl")
include("./profiles.jl")
include("./cib.jl")
include("./co_broadband.jl")
include("./lrg.jl")
include("./radio.jl")

export get_cosmology, read_halo_catalog_hdf5, sort_halo_catalog
export Radio_Sehgal2009, CIB_Planck2013, CIB_Scarfy, CO_CROWNED, LRG_Yuan23
export paint!, generate_sources, process_sources, profile_grid, profile_paint!, profileworkspace
export build_interpolator, Battaglia16ThermalSZProfile

end # module
