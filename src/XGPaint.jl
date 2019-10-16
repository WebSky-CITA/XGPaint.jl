module XGPaint


include("model.jl")
include("cib.jl")


# function fill_halovars!(
#         x, y, z, # inputs
#         redshift_result, dist_result) # outputs
#     """
#     This function computes distance and redshift in parallel.
#     """
#
#     N_halos = size(x,1)
#     @sync @distributed for i = 1:N_halos
#         dist_result[i] = sqrt(x[i]^2 + y[i]^2 + z[i]^2)
#         redshift_result[i] = r2z(dist_result[i])
#     end
# end
#
# function hod_shang(cen_result, sat_result, sat_bar_result,
#         halo_mass::SharedArray, par)
#     # computes shang HOD and generates a Poisson draw
#     min_lm = log(minimum(halo_mass))
#     max_lm = log(maximum(halo_mass))
#     logM_to_N = build_shang_interp(min_lm, max_lm, par)
#     N_halos = size(halo_mass,1)
#
#     @sync @distributed for i = 1:N_halos
#         cen_result[i] = 1
#
#         sat_bar_result[i] = logM_to_N(log(halo_mass[i]))
#         sat_result[i] = pois_rand(convert(Float64, sat_bar_result[i]))
#     end
#
# end


end # module
