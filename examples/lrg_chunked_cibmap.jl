# hackish script for generating LRG/CIB maps after the fact
using XGPaint
using Healpix
using JLD2, FileIO, CodecZlib
using HDF5
using DelimitedFiles

cosmo = get_cosmology(h=0.677f0, OmegaM=0.310f0)
model = LRG_Yuan22{Float32}(nside=4096)

freqs = ["353", "545", "857",
    "18.7", "21.6", "24.5", "27.3", "30.0", "35.9", "41.7", "44.0", "47.4",
    "63.9", "67.8", "70.0", "73.7", "79.6", "90.2", "100", "111", "129", "143",
    "153", "164", "189", "210", "217", "232", "256", "275", "294", "306", "314",
    "340", "375", "409", "467", "525", "584", "643", "729", "817",
    "906", "994", "1080"
]

function write_chunk(
                     lrg_dir, cib_dir, chunk_index, model, cosmo, freqs)
    m = HealpixMap{Float64,RingOrder}(model.nside)
    jldopen(joinpath(lrg_dir, "sources/cen_chunk$(chunk_index).jld2"), "r") do file
        println("loading centrals")
        lrg_cen = file["LRG"]
        println(sum(lrg_cen),"/",length(lrg_cen))
        hp_ind_lrgcen = file["hp_ind_cen"][lrg_cen]
        # process centrals for this frequency
        for freq in freqs
            print(freq,"...")
            fill!(m.pixels, zero(Float32))
            h5open(joinpath(cib_dir, "sources/cen_chunk$(chunk_index)_flux_$(freq).h5"), "r") do file2
                fluxes_lrgcen = read(file2["flux"])[lrg_cen]
                Threads.@threads for i in 1:length(hp_ind_lrgcen)
                    m.pixels[hp_ind_lrgcen[i]] += fluxes_lrgcen[i]
                end
            end
            # save maps
            filename = joinpath(lrg_dir, "lrg_cen_cib_$(freq).fits")

            if chunk_index > 1
                m0 = Healpix.readMapFromFITS(filename, 1, Float32)
                m.pixels = m.pixels + m0.pixels
            end
            Healpix.saveToFITS(m, "!$(filename)", typechar="D")
        end
    end
    println("chunk done")
end

function run_all_chunks(lrg_dir, cib_dir, freqs, N_chunks)
    for chunk_index in 1:N_chunks
        println("Chunk ", chunk_index, "/", N_chunks)
        write_chunk(lrg_dir, cib_dir, chunk_index, model, cosmo, freqs)
    end
end
## compute on all chunks, on all halos
scratch_dir = "/home/dongwooc/scratchspace/lrg_sources"
scratch2dir = "/home/dongwooc/scratchspace/cib_co_sources"
run_all_chunks(scratch_dir, scratch2dir, freqs, 4)
