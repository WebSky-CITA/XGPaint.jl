
using XGPaint

using Cosmology
using Healpix
import Future
import Distributions
using Random

cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = CIBModel{Float64}()
res = Resolution(model.nside)
##
halo_pos_serial, halo_mass_serial = read_halo_catalog_hdf5(
    "/home/zequnl/websky_halos-light.hdf5")
##
function generate_sources_threaded(
        # model parameters
        model::CIBModel, cosmo::Cosmology.FlatLCDM{T}, Healpix_res::Resolution,
        # halo arrays
        halo_pos::Array{T,2}, halo_mass::Array{T,1};
        verbose=true) where T

    N_halos = size(halo_mass, 1)
    interp = get_interpolators( model, cosmo,
        minimum(halo_mass), maximum(halo_mass))

    verbose && println("Allocating arrays.")

    hp_ind_cen = Array{Int64}(undef, N_halos)  # healpix index of halo
    lum_cen = Array{T}(undef, N_halos)  # Lum of central w/o ν-dependence
    redshift_cen = Array{T}(undef, N_halos)
    dist_cen = Array{T}(undef, N_halos)
    n_sat_bar = Array{T}(undef, N_halos)
    n_sat_bar_result = Array{Int32}(undef, N_halos)

    # STEP 1: compute central properties -----------------------------------
    verbose && println("Processing centrals on $(Threads.nthreads()) threads.")
    Threads.@threads for i = 1:N_halos
        # location information for centrals
        hp_ind_cen[i] = Healpix.vec2pixRing(Healpix_res,
            halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
        dist_cen[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
        redshift_cen[i] = interp.r2z(dist_cen[i])

        # compute HOD
        n_sat_bar[i] = interp.hod_shang(log(halo_mass[i]))
        n_sat_bar_result[i] = rand(Distributions.PoissonADSampler(
            Float64(n_sat_bar[i])))

        # get central luminosity
        lum_cen[i] = sigma_cen(halo_mass[i], model) * (
            one(T) + redshift_cen[i])^model.shang_eta
    end

    # STEP 2: Generate satellite arrays -----------------------------
    # set up indices for satellites
    cumsat = cumsum(n_sat_bar_result)
    prepend!(cumsat, 0)
    total_n_sat = cumsat[end]
    # result arrays
    hp_ind_sat = Array{Int64}(undef, total_n_sat)  # healpix index of halo
    lum_sat = Array{T}(undef, total_n_sat)  # Lum of central w/o ν-dependence
    redshift_sat = Array{T}(undef, total_n_sat)

    verbose && println("Processing $(total_n_sat) satellites.")
    Threads.@threads for i = 1:N_halos
        r_cen = m2r(halo_mass[i], cosmo)
        c_cen = mz2c(halo_mass[i], redshift_cen[i], cosmo)
        for j in 1:n_sat_bar_result[i]
            # i is central halo index, j is index of satellite within each halo
            i_sat = cumsat[i]+j # index of satellite in satellite arrays

            log_msat_inner = max(log(rand(T)), T(-7.0))
            r_sat = r_cen * interp.c_lnm2r(
                c_cen, log_msat_inner) * T(200.0^(-1/3))
            m_sat =  interp.muofn(rand(T) * n_sat_bar[i]) * halo_mass[i]

            phi = random_phi(Float32)
            theta = random_theta(Float32)
            x_sat = halo_pos[1,i] + r_sat * sin(theta) * cos(phi)
            y_sat = halo_pos[2,i] + r_sat * sin(theta) * sin(phi)
            z_sat = halo_pos[3,i] + r_sat * cos(theta)
            d_sat = sqrt(x_sat^2 + y_sat^2 + z_sat^2)
            redshift_sat[i_sat] = interp.r2z(d_sat)

            lum_sat[i_sat] = sigma_cen(m_sat, model) * (
                1+redshift_sat[i_sat])^model.shang_eta
            hp_ind_sat[i_sat] = Healpix.vec2pixRing(
                Healpix_res, x_sat, y_sat, z_sat)
        end
    end
    
    return hp_ind_cen, lum_cen, redshift_cen, hp_ind_sat, lum_sat, redshift_sat
end


##
@time generate_sources_threaded(model, cosmo, res, halo_pos_serial, halo_mass_serial)
