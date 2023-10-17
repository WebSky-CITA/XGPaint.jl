using XGPaint
using SparseArrays
using Healpix
using HDF5
using JLD2, FileIO, CodecZlib
using DelimitedFiles
using Random
Random.seed!(3)

halo_pos, halo_mass = read_halo_catalog_hdf5("/fs/lustre/cita/zack/projects/websky/websky_halos-light.hdf5")
cosmo = get_cosmology(h=0.677f0, OmegaM=0.310f0)
modelCIB = CIB_Scarfy{Float32}()
model_CO = CO_CROWNED{Float32}()
freqs = [
    "18.7", "21.6", "24.5", "27.3", "30.0", "35.9", "41.7", "44.0", "47.4",
    "63.9", "67.8", "70.0", "73.7", "79.6", "90.2", "100", "111", "129", "143",
    "153", "164", "189", "210", "217", "232", "256", "275", "294", "306", "314",
    "340", "353", "375", "409", "467", "525", "545", "584", "643", "729", "817",
    "857", "906", "994", "1080"
]

bandpass_PA6_090 = open(readdlm,"/home/dongwooc/lustrespace/tilec/data/PA6_avg_passband_f090_wErr_trunc_20200505.txt")
bandpass_PA6_150 = open(readdlm,"/home/dongwooc/lustrespace/tilec/data/PA6_avg_passband_f150_wErr_trunc_20200505.txt")
bandpass_PA4_220 = open(readdlm,"/home/dongwooc/lustrespace/tilec/data/PA4_avg_passband_f220_wErr_trunc_20200505.txt")

function write_chunk(
                     output_dir, chunk_index, modelCIB, model_CO, cosmo,
                     pos, mass, freqs)
    sourcesCIB = generate_sources(modelCIB, cosmo, pos, mass);
    sources_CO = process_sources(model_CO, sourcesCIB, modelCIB);
    fluxes_cen = Array{Float32,1}(undef,sourcesCIB.N_cen)
    fluxes_sat = Array{Float32,1}(undef,sourcesCIB.N_sat)
    println("writing info for ",sourcesCIB.N_cen," centrals")
    h5open(joinpath(output_dir, "sources/cen_chunk$(chunk_index).h5"), "w") do file
        write(file, "redshift", sourcesCIB.redshift_cen)
        write(file, "theta", sourcesCIB.theta_cen)
        write(file, "phi", sourcesCIB.phi_cen)
        write(file, "LIR", sources_CO.LIR_cen)
        write(file, "LCO", sources_CO.LcoJ_cen)
        write(file, "LCI", sources_CO.LCI_cen)
    end
    println("writing info for ",sourcesCIB.N_sat," satellites")
    h5open(joinpath(output_dir, "sources/sat_chunk$(chunk_index).h5"), "w") do file
        write(file, "redshift", sourcesCIB.redshift_sat)
        write(file, "theta", sourcesCIB.theta_sat)
        write(file, "phi", sourcesCIB.phi_sat)
        write(file, "LIR", sources_CO.LIR_sat)
        write(file, "LCO", sources_CO.LcoJ_sat)
        write(file, "LCI", sources_CO.LCI_sat)
    end
    m = HealpixMap{Float64,RingOrder}(modelCIB.nside)
    mCO090 = HealpixMap{Float64,RingOrder}(model_CO.nside)
    mCO150 = HealpixMap{Float64,RingOrder}(model_CO.nside)
    mCO220 = HealpixMap{Float64,RingOrder}(model_CO.nside)
    mCI090 = HealpixMap{Float64,RingOrder}(model_CO.nside)
    mCI150 = HealpixMap{Float64,RingOrder}(model_CO.nside)
    mCI220 = HealpixMap{Float64,RingOrder}(model_CO.nside)
    println("***  90 GHz ***")
    fluxCO_090_cen, fluxCI_090_cen, fluxCO_090_sat, fluxCI_090_sat = XGPaint.paint!(mCO090, mCI090, bandpass_PA6_090, model_CO, sources_CO)
    println("*** 150 GHz ***")
    fluxCO_150_cen, fluxCI_150_cen, fluxCO_150_sat, fluxCI_150_sat = XGPaint.paint!(mCO150, mCI150, bandpass_PA6_150, model_CO, sources_CO)
    println("*** 220 GHz ***")
    fluxCO_220_cen, fluxCI_220_cen, fluxCO_220_sat, fluxCI_220_sat = XGPaint.paint!(mCO220, mCI220, bandpass_PA4_220, model_CO, sources_CO)
    println("writing CO maps ...")
    filename = joinpath(output_dir, "co_090.fits")
    if chunk_index > 1
        m0 = Healpix.readMapFromFITS(filename, 1, Float32)
        mCO090.pixels = mCO090.pixels + m0.pixels
    end
    Healpix.saveToFITS(mCO090, "!$(filename)", typechar="D")
    filename = joinpath(output_dir, "co_150.fits")
    if chunk_index > 1
        m0 = Healpix.readMapFromFITS(filename, 1, Float32)
        mCO150.pixels = mCO150.pixels + m0.pixels
    end
    Healpix.saveToFITS(mCO150, "!$(filename)", typechar="D")
    filename = joinpath(output_dir, "co_220.fits")
    if chunk_index > 1
        m0 = Healpix.readMapFromFITS(filename, 1, Float32)
        mCO220.pixels = mCO220.pixels + m0.pixels
    end
    Healpix.saveToFITS(mCO220, "!$(filename)", typechar="D")
    println("writing CI maps ...")
    filename = joinpath(output_dir, "ci_090.fits")
    if chunk_index > 1
        m0 = Healpix.readMapFromFITS(filename, 1, Float32)
        mCI090.pixels = mCI090.pixels + m0.pixels
    end
    Healpix.saveToFITS(mCI090, "!$(filename)", typechar="D")
    filename = joinpath(output_dir, "ci_150.fits")
    if chunk_index > 1
        m0 = Healpix.readMapFromFITS(filename, 1, Float32)
        mCI150.pixels = mCI150.pixels + m0.pixels
    end
    Healpix.saveToFITS(mCI150, "!$(filename)", typechar="D")
    filename = joinpath(output_dir, "ci_220.fits")
    if chunk_index > 1
        m0 = Healpix.readMapFromFITS(filename, 1, Float32)
        mCI220.pixels = mCI220.pixels + m0.pixels
    end
    Healpix.saveToFITS(mCI220, "!$(filename)", typechar="D")
    println("writing CO fluxes ...")
    for J in 1:7
        jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCO$(J)_090.jld2"), "w") do file
            write(file, "flux", sparse(fluxCO_090_cen[:,J]))
        end
        jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCO$(J)_150.jld2"), "w") do file
            write(file, "flux", sparse(fluxCO_150_cen[:,J]))
        end
        jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCO$(J)_220.jld2"), "w") do file
            write(file, "flux", sparse(fluxCO_220_cen[:,J]))
        end
        jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCO$(J)_090.jld2"), "w") do file
            write(file, "flux", sparse(fluxCO_090_sat[:,J]))
        end
        jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCO$(J)_150.jld2"), "w") do file
            write(file, "flux", sparse(fluxCO_150_sat[:,J]))
        end
        jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCO$(J)_220.jld2"), "w") do file
            write(file, "flux", sparse(fluxCO_220_sat[:,J]))
        end
    end
    println("writing CI fluxes ...")
    jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCI_090.jld2"), "w") do file
        write(file, "flux", sparse(fluxCI_090_cen))
    end
    jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCI_150.jld2"), "w") do file
        write(file, "flux", sparse(fluxCI_150_cen))
    end
    jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCI_220.jld2"), "w") do file
        write(file, "flux", sparse(fluxCI_220_cen))
    end
    jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCI_090.jld2"), "w") do file
        write(file, "flux", sparse(fluxCI_090_sat))
    end
    jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCI_150.jld2"), "w") do file
        write(file, "flux", sparse(fluxCI_150_sat))
    end
    jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCI_220.jld2"), "w") do file
        write(file, "flux", sparse(fluxCI_220_sat))
    end
    for freq in freqs

        XGPaint.paint!(m, parse(Float32, freq) * 1.0f9, modelCIB, sourcesCIB,
            fluxes_cen, fluxes_sat)
        # save fluxes
        h5open(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_flux_$(freq).h5"), "w") do file
            write(file, "flux", fluxes_cen)
        end
        h5open(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_flux_$(freq).h5"), "w") do file
            write(file, "flux", fluxes_sat)
        end
        # save maps
        filename = joinpath(output_dir, "cib_$(freq).fits")

        if chunk_index > 1
            m0 = Healpix.readMapFromFITS(filename, 1, Float32)
            m.pixels = m.pixels + m0.pixels
        end
        Healpix.saveToFITS(m, "!$(filename)", typechar="D")
    end
end

function run_all_chunks(output_dir, halo_pos, halo_mass, freqs; N_chunks=4)
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
        write_chunk(output_dir, chunk_index, modelCIB, model_CO, cosmo,
            pos, mass, freqs)
    end
    # save maps
    println("writing CO+CI maps ...")
    m = HealpixMap{Float64,RingOrder}(model_CO.nside)
    filename = joinpath(output_dir, "co_ci_090.fits")
    m0 = Healpix.readMapFromFITS(joinpath(output_dir, "co_090.fits"),
                                  1, Float32)
    m1 = Healpix.readMapFromFITS(joinpath(output_dir, "ci_090.fits"),
                                  1, Float32)
    m.pixels = m0.pixels + m1.pixels
    Healpix.saveToFITS(m, "!$(filename)", typechar="D")
    filename = joinpath(output_dir, "co_ci_150.fits")
    m0 = Healpix.readMapFromFITS(joinpath(output_dir, "co_150.fits"),
                                  1, Float32)
    m1 = Healpix.readMapFromFITS(joinpath(output_dir, "ci_150.fits"),
                                  1, Float32)
    m.pixels = m0.pixels + m1.pixels
    Healpix.saveToFITS(m, "!$(filename)", typechar="D")
    filename = joinpath(output_dir, "co_ci_220.fits")
    m0 = Healpix.readMapFromFITS(joinpath(output_dir, "co_220.fits"),
                                  1, Float32)
    m1 = Healpix.readMapFromFITS(joinpath(output_dir, "ci_220.fits"),
                                  1, Float32)
    m.pixels = m0.pixels + m1.pixels
    Healpix.saveToFITS(m, "!$(filename)", typechar="D")
end
## compute on all chunks, on all halos
scratch_dir = "/home/dongwooc/projectscratchspace/cib_co_sources_scarfy_refactor"
println("SCRATCH: ", scratch_dir)
mkpath(scratch_dir)
mkpath(joinpath(scratch_dir, "sources"))
run_all_chunks(scratch_dir, halo_pos, halo_mass, freqs)
