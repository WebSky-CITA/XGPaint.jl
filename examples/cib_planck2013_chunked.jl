using XGPaint
using Healpix

## Load halos from HDF5 files, establish a CIB model and cosmology
halo_pos, halo_mass = read_halo_catalog_hdf5(
    "/global/cscratch1/sd/xzackli/websky_halos-light.hdf5")

cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = CIB_Planck2013{Float32}(nside=8192)

## Write one chunk to disk
function write_chunk(output_dir, chunk_index, model, cosmo, halo_pos, halo_mass, freqs)
    # Allocate some arrays and fill them up for centrals and satellites
    @time sources = generate_sources(model, cosmo, halo_pos, halo_mass);

    # Deposit the sources into maps
    fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
    fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)
    m = Map{Float64,RingOrder}(model.nside)

    # loop over all frequencies and paint sources to appropriate freq map
    @time begin
        for freq in freqs
            XGPaint.paint!(m, parse(Float32, freq) * 1.0f9, model, sources,
                fluxes_cen, fluxes_sat)

            # read from disk if not the first chunk
            filename = "$(output_dir)/cib_$(freq).fits"
            if chunk_index > 1
                m0 = Healpix.readMapFromFITS(filename, 1, Float64)
                m.pixels = m.pixels + m0
            end
            Healpix.saveToFITS(m, "!$(filename)")
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
        println("Chunk ", chunk_index, "/", length(chunks))
        left_ind, right_ind = chunks[chunk_index]

        pos = @view halo_pos[:, left_ind:right_ind]
        mass = @view halo_mass[left_ind:right_ind]
        write_chunk(output_dir, chunk_index, model, cosmo,
            pos, mass, freqs)
    end
end
## compute on all chunks, on all halos

# freqs = ["143"]
freqs = [
    "18.7", "21.6", "24.5", "27.3", "30.0", "35.9", "41.7", "44.0", "47.4",
    "63.9", "67.8", "70.0", "73.7", "79.6", "90.2", "100", "111", "129", "143",
    "153", "164", "189", "210", "217", "232", "256", "275", "294", "306", "314",
    "340", "353", "375", "409", "467", "525", "545", "584", "643", "729", "817",
    "857", "906", "994", "1080"
]
scratch_dir = "/global/cscratch1/sd/xzackli/cib/"
run_all_chunks(scratch_dir, halo_pos, halo_mass, freqs; N_chunks=2)

##
