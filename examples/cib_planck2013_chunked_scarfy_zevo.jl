using XGPaint
using Healpix
using HDF5

## Load halos from HDF5 files, establish a CIB model and cosmology
@time halo_pos, halo_mass = read_halo_catalog_hdf5("/fs/lustre/cita/zack/projects/websky/websky_halos-light.hdf5")
cosmo = get_cosmology(h=0.677f0, OmegaM=0.310f0)
model = CIB_Planck2013{Float32}(nside=4096,z_evo="scarfy",quench=true)

## Write one chunk to disk

function write_chunk(
    output_dir, chunk_index, model, cosmo,
    pos, mass, freqs)
    # Allocate some arrays and fill them up for centrals and satellites
    @time sources = generate_sources(model, cosmo, pos, mass);

    # Deposit the sources into maps
    fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
    fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)
    m = HealpixMap{Float64,RingOrder}(model.nside)

    #h5open(joinpath(output_dir, "sources/cen_chunk$(chunk_index).h5"), "w") do file
    #    write(file, "redshift", sources.redshift_cen)
    #    write(file, "theta", sources.theta_cen)
    #    write(file, "phi", sources.phi_cen)
    #end
    #h5open(joinpath(output_dir, "sources/sat_chunk$(chunk_index).h5"), "w") do file
    #    write(file, "redshift", sources.redshift_sat)
    #    write(file, "theta", sources.theta_sat)
    #    write(file, "phi", sources.phi_sat)
    #end

    # loop over all frequencies and paint sources to appropriate freq map
    @time begin
        for freq in freqs

            XGPaint.paint!(m, parse(Float32, freq) * 1.0f9, model, sources,
                fluxes_cen, fluxes_sat)

            # save sources with mass, redshift, angles
            #h5open(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_flux_$(freq).h5"), "w") do file
            #    write(file, "flux", fluxes_cen)
            #end
            #h5open(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_flux_$(freq).h5"), "w") do file
            #    write(file, "flux", fluxes_sat)
            #end

            filename = joinpath(output_dir, "cib_$(freq).fits")

            if chunk_index > 1
                m0 = Healpix.readMapFromFITS(filename, 1, Float32)
                m.pixels = m.pixels + m0.pixels
            end
            Healpix.saveToFITS(m, "!$(filename)", typechar="D")
        end
    end

end

##
function run_all_chunks(output_dir, halo_pos, halo_mass, freqs; N_chunks=2)
    # provide views into halo positions and masses for chunks of the halos
    N_halos = size(halo_mass, 1)
    chunksize = trunc(Integer, N_halos / N_chunks + 1)
    chunks = chunk(N_halos, chunksize)
    for chunk_index in 1:length(chunks)
    left_ind, right_ind = chunks[chunk_index]
        println("Chunk ", chunk_index, "/", length(chunks),
            " ", left_ind, " ", right_ind)
        pos = @view halo_pos[:, left_ind:right_ind]
        mass = @view halo_mass[left_ind:right_ind]
        write_chunk(output_dir, chunk_index, model, cosmo,
            pos, mass, freqs)
    end
end
## compute on all chunks, on all halos

"""
freqs = [
    "18.7", "21.6", "24.5", "27.3", "30.0", "35.9", "41.7", "44.0", "47.4",
    "63.9", "67.8", "70.0", "73.7", "79.6", "90.2", "100", "111", "129", "143",
    "153", "164", "189", "210", "217", "232", "256", "275", "294", "306", "314",
    "340", "353", "375", "409", "467", "525", "545", "584", "643", "729", "817",
    "857", "906", "994", "1080"
]
"""
freqs = ["143", "857", "217", "353", "545"]

scratch_dir = "/home/dongwooc/scratchspace/cib_scarfy_test"
println("SCRATCH: ", scratch_dir)
mkpath(scratch_dir)
mkpath(joinpath(scratch_dir, "sources"))
run_all_chunks(scratch_dir, halo_pos, halo_mass, freqs)

##
