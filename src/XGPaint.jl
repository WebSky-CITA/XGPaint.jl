
module XGPaint

using Distributed

include("./model.jl")
include("./cib.jl")
include("./util.jl")


"""
Compute sources for CIB.
"""
function generate_sources(
        # model parameters
        model::CIBModel, cosmo::Cosmology.FlatLCDM{T}, interp::NamedTuple,
        Healpix_res::Resolution,
        # halo arrays
        halo_pos::SharedArray{T,2}, halo_mass::SharedArray{T,1}) where T

    N_halos = size(halo_mass, 1)

    # result arrays
    hp_ind_cen = SharedArray{Int64}(N_halos)  # healpix index of halo
    lum_cen = SharedArray{T}(N_halos)  # Lum of central w/o ν-dependence
    redshift_cen = SharedArray{T}(N_halos)
    dist_cen = SharedArray{T}(N_halos)
    sat_bar = SharedArray{T}(N_halos)
    sat_bar_result = SharedArray{Int32}(N_halos)

    # STEP 1: compute central properties -----------------------------------
    @sync @distributed for i = 1:N_halos
        # location information for centrals
        hp_ind_cen[i] = Healpix.vec2pixRing(Healpix_res,
            halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
        dist_cen[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
        redshift_cen[i] = interp.r2z(dist_cen[i])

        # compute HOD
        sat_bar[i] = interp.hod_shang(log(halo_mass[i]))
        sat_bar_result[i] = pois_rand(convert(Float64, sat_bar[i]))

        # get central luminosity
        lum_cen[i] = sigma_cen(halo_mass[i], model) * (
            one(T) + redshift_cen[i])^model.shang_eta
    end

    return sat_bar_result
end



end # module
