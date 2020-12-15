using XGPaint
using Healpix
using HDF5
import ThreadsX

ENV["SCRATCH"] = "/home/zequnl/scratch/"
@time halo_pos, halo_mass = read_halo_catalog_hdf5(
    joinpath(ENV["SCRATCH"],"websky_halos-light.hdf5"));

## Load halos from HDF5 files, establish a CIB model and cosmology
cosmo = get_cosmology(Float32, h = 0.7, OmegaM = 0.25)
radio_model = Radio_Sehgal2009{Float32}(a_0 = -1)
@time begin
    sources = generate_sources(radio_model, cosmo, halo_pos, halo_mass);
end;

##
# freqs = [
#     "18.7", "21.6", "24.5", "27.3", "30.0", "35.9", "41.7", "44.0", "47.4", "63.9", "67.8",
#     "70.0", "73.7", "79.6", "90.2", "100", "111", "129", "143", "153", "164", "189", "210", 
#     "217", "232", "256", "275", "294", "306", "314", "340", "353", "375", "409", "467", 
#     "525", "545", "584", "643", "729", "817", "857", "906", "994", "1080"
# ]

freqs = ["143.0"]

m = Map{Float64,RingOrder}(radio_model.nside)

function generate_maps()
    scratch_dir = ENV["SCRATCH"]
    println("SCRATCH: ", scratch_dir)

    for freq in freqs
        @time begin
            flux_I, redshift_I, θ_I, ϕ_I, flux_II, redshift_II, θ_II, ϕ_II = paint!(
                m, parse(Float32, freq) * 1.0f9, radio_model, sources,
                return_fluxes=true)

            flux = vcat(flux_I, flux_II)
            redshift = vcat(redshift_I, redshift_II)
            theta = vcat(θ_I, θ_II)
            phi = vcat(ϕ_I, ϕ_II)
            # perm = sortperm(flux, rev=true, alg=ThreadsX.MergeSort)

            h5open(joinpath(scratch_dir, "catalog_$(freq).h5"), "w") do file
                write(file, "flux", flux)
                write(file, "redshift", redshift)
                write(file, "theta", theta)
                write(file, "phi", phi)
            end
            # map_path = joinpath(scratch_dir, "radio$(freq).fits")
            # Healpix.saveToFITS(m, "!$(map_path)")
        end
    end
end

generate_maps()
##
