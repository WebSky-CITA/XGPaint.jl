##
using Distributed
addprocs(16)
##

@everywhere begin
    using XGPaint

    using SharedArrays
    using Cosmology
    using Healpix

    cosmo = XGPaint.get_cosmology(h=0.7f0, OmegaM=0.25f0)
    model = XGPaint.CIBModel{Float32}()
    interp_CIB = XGPaint.get_interpolators( model, cosmo, 1.0f12, 3.0f15 )
    res = Resolution(model.nside)
end

##
# only the master has these
#
halo_pos_serial, halo_mass_serial = read_halo_catalog_hdf5(
    "/home/zequnl/websky_halos-light.hdf5")

# halo_pos_serial = Float32.(rand(3,100))
# halo_mass_serial = Float32.(1e13 .* (1 .+ rand(100)))
pos = convert(SharedArray, halo_pos_serial)
mass = convert(SharedArray, halo_mass_serial)
halo_mass_serial, halo_pos_serial = 0, 0
@everywhere Base.GC.gc()

##
XGPaint.generate_sources(model, cosmo, interp_CIB, res, pos, mass)
##
# generate_sources(model, cosmo, interp_CIB, res, pos, mass)

using PoissonRandom
using XGPaint

using SharedArrays
using Cosmology
using Healpix
import Future
import Distributions
using Random

cosmo = XGPaint.get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = XGPaint.CIBModel{Float32}()
interp_CIB = XGPaint.get_interpolators( model, cosmo, 1.0f12, 3.0f15 )
res = Resolution(model.nside)
##
halo_pos_serial, halo_mass_serial = read_halo_catalog_hdf5(
    "/home/zequnl/websky_halos-light.hdf5")
##
function generate_sources_threaded(
        # model parameters
        model::CIBModel, cosmo::Cosmology.FlatLCDM{T}, interp::NamedTuple,
        Healpix_res::Resolution,
        # halo arrays
        halo_pos::Array{T,2}, halo_mass::Array{T,1}) where T

    N_halos = size(halo_mass, 1)

    # result arrays
    hp_ind_cen = Array{Int64}(undef, N_halos)  # healpix index of halo
    lum_cen = Array{T}(undef, N_halos)  # Lum of central w/o ν-dependence
    redshift_cen = Array{T}(undef, N_halos)
    dist_cen = Array{T}(undef, N_halos)
    sat_bar = Array{T}(undef, N_halos)
    sat_bar_result = Array{Int32}(undef, N_halos)

    # STEP 1: compute central properties -----------------------------------
    nchunks = Threads.nthreads() * 10
    chunksize = Int64(ceil(N_halos / nchunks))
    chunks = [(chunksize*(nc-1)+1, min(chunksize*nc, N_halos) ) for nc in 1:nchunks]
    println(chunks)

    Threads.@threads for thread_i = 1:nchunks
        li, ri = chunks[thread_i]
        for i in li:ri
            # location information for centrals
            hp_ind_cen[i] = Healpix.vec2pixRing(Healpix_res,
                halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
            dist_cen[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
            redshift_cen[i] = interp.r2z(dist_cen[i])

            # compute HOD
            sat_bar[i] = interp.hod_shang(log(halo_mass[i]))
            sat_bar_result[i] = rand(Distributions.PoissonADSampler(Float64(sat_bar[i])))

            # get central luminosity
            lum_cen[i] = sigma_cen(halo_mass[i], model) * (
                one(T) + redshift_cen[i])^model.shang_eta
        end
    end

    return sat_bar
end


##
generate_sources_threaded(model, cosmo, interp_CIB, res, halo_pos_serial, halo_mass_serial)
