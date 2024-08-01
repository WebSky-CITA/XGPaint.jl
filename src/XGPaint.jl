
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
using DelimitedFiles

import PhysicalConstants.CODATA2018 as constants
const M_sun = 1.98847e30u"kg"
const T_cmb =  2.725 * u"K"
const P_e_factor = constants.σ_e / (constants.m_e * constants.c_0^2)

include("./util.jl")
include("./model.jl")
include("./profiles.jl")
include("./profiles_rsz.jl")
include("./profiles_szp.jl")
include("./tau.jl")
include("./cib.jl")
include("./co_broadband.jl")
include("./lrg.jl")
include("./radio.jl")

export get_cosmology, read_halo_catalog_hdf5, sort_halo_catalog
export Radio_Sehgal2009, CIB_Planck2013, CIB_Scarfy, CO_CROWNED, LRG_Yuan23
export paint!, generate_sources, process_sources, profile_grid, profile_paint!, profileworkspace, paint_szp!, profile_grid_szp, profile_paint_szp!, paint_rsz!, profile_grid_rsz, profile_paint_rsz!
export build_interpolator, Battaglia16ThermalSZProfile, Battaglia16RelativisticSZProfile, build_interpolator_szp, build_interpolator_rsz
export Battaglia16SZPackProfile, nu_to_X, X_to_nu, BattagliaTauProfile

end # module
