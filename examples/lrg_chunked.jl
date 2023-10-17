using XGPaint
using Healpix
using JLD2, FileIO, CodecZlib
using HDF5
using DelimitedFiles
using Random
Random.seed!(3)

hdata = h5open("/home/dongwooc/projectscratchspace/websky_halos_rewrite/websky_halos-lesslight_20230612.h5","r") do file
   (x = read(file,"x"), y = read(file,"y"), z = read(file,"z"), m = read(file,"M200m"), vrad = read(file,"vrad") )
end
halo_pos = vcat(hdata.x',hdata.y',hdata.z')
halo_mass = hdata.m
halo_vrad = hdata.vrad
cosmo = get_cosmology(h=0.677f0, OmegaM=0.310f0)
model = LRG_Yuan23{Float32}(nside=4096)

lrg_redshifts = [
    "0.35", "0.40", "0.45", "0.50", "0.55",
    "0.60", "0.65", "0.70", "0.75", "0.80"
]

function write_chunk(
                     output_dir, chunk_index, model, cosmo,
                     pos, mass, vrad, lrg_redshifts)
    sources = generate_sources(model, cosmo, pos, mass, vrad);
    jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index).jld2"), "w") do file
        file["redshift"]=sources.redshift_cen
        file["theta"]=sources.theta_cen
        file["phi"]=sources.phi_cen
        write(file,"LRG",sources.lrg_cen)
        file["m200c"]=sources.m200c_cen
        file["hp_ind_cen"]=sources.hp_ind_cen
    end
    jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index).jld2"), "w") do file
        file["redshift"]=sources.redshift_sat
        file["theta"]=sources.theta_sat
        file["phi"]=sources.phi_sat
        write(file,"LRG",sources.lrg_sat)
        file["m200m"]=sources.m200m_sat
        file["m200c"]=sources.m200c_sat
        file["parent"]=sources.parent_sat
        file["hp_ind_sat"]=sources.hp_ind_sat
    end
    m = HealpixMap{Float64,RingOrder}(model.nside)
    for zlo in lrg_redshifts
        XGPaint.paint!(m, model, sources,
            parse(Float32, zlo), parse(Float32, zlo)+0.05f0)
        # save maps
        filename = joinpath(output_dir, "lrg_$(zlo).fits")

        if chunk_index > 1
            m0 = Healpix.readMapFromFITS(filename, 1, Float32)
            m.pixels = m.pixels + m0.pixels
        end
        Healpix.saveToFITS(m, "!$(filename)", typechar="D")
    end
end

## the rest of this is more or less the same as cib_planck2013_chunked.jl
function run_all_chunks(output_dir, halo_pos, halo_mass, halo_vrad, lrg_redshifts; N_chunks=4)
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
        vrad = @view halo_vrad[left_ind:right_ind]
        write_chunk(output_dir, chunk_index, model, cosmo,
            pos, mass, vrad, lrg_redshifts)
    end
end
## compute on all chunks, on all halos
scratch_dir = "/home/dongwooc/projectscratchspace/lrg_sources_vrad_m200c"
println("SCRATCH: ", scratch_dir)
mkpath(scratch_dir)
mkpath(joinpath(scratch_dir, "sources"))
run_all_chunks(scratch_dir, halo_pos, halo_mass, halo_vrad, lrg_redshifts)
