# OBSOLETE! DO NOT USE WITH CURRENT XGPAINT
# hackish script for generating LRG maps after the fact
using XGPaint
using Healpix
using JLD2, FileIO, CodecZlib
using DelimitedFiles

cosmo = get_cosmology(h=0.677f0, OmegaM=0.310f0)
model = LRG_Yuan22{Float32}(nside=4096)

zlo = 0.64
zhi = 0.86

function write_chunk(
                     output_dir, chunk_index, model, cosmo, zlo, zhi)
    m = HealpixMap{Float64,RingOrder}(model.nside)
    fill!(m.pixels, zero(Float32))
    jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index).jld2"), "r") do file
        redshift_cen = file["redshift"]
        lrg_cen = file["LRG"]
        hp_ind_cen = file["hp_ind_cen"]
        # process centrals for this frequency
        Threads.@threads for i in 1:length(redshift_cen)
            if (lrg_cen[i]) && (redshift_cen[i] < zhi) && (redshift_cen[i] > zlo)
                m.pixels[hp_ind_cen[i]] += 1
            end
        end
    end
    jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index).jld2"), "r") do file
        redshift_sat = file["redshift"]
        lrg_sat = file["LRG"]
        hp_ind_sat = file["hp_ind_sat"]
        # process satellites for this frequency
        Threads.@threads for i in 1:length(redshift_sat)
            if (lrg_sat[i]) && (redshift_sat[i] < zhi) && (redshift_sat[i] > zlo)
                m.pixels[hp_ind_sat[i]] += 1
            end
        end
    end
    # save maps
    filename = joinpath(output_dir, "lrg_$(zlo)_$(zhi).fits")

    if chunk_index > 1
        m0 = Healpix.readMapFromFITS(filename, 1, Float32)
        m.pixels = m.pixels + m0.pixels
    end
    Healpix.saveToFITS(m, "!$(filename)", typechar="D")
end

function run_all_chunks(output_dir, zlo, zhi, N_chunks)
    for chunk_index in 1:N_chunks
        println("Chunk ", chunk_index, "/", N_chunks)
        write_chunk(output_dir, chunk_index, model, cosmo, zlo, zhi)
    end
end
## compute on all chunks, on all halos
scratch_dir = "/home/dongwooc/scratchspace/lrg_sources"
println("SCRATCH: ", scratch_dir)
mkpath(scratch_dir)
mkpath(joinpath(scratch_dir, "sources"))
run_all_chunks(scratch_dir, zlo, zhi, 4)
