using XGPaint
using Healpix
using HDF5
using JLD2, FileIO, CodecZlib
using DelimitedFiles
using Random
Random.seed!(3)

halo_pos, halo_mass = read_halo_catalog_hdf5("/fs/lustre/cita/zack/projects/websky/websky_halos-light.hdf5")
cosmo = get_cosmology(h=0.677f0, OmegaM=0.310f0)
model = CIB_Planck2013{Float32}(z_evo="scarfy",m_evo="scarfy")
println(model.z_evo)
println(model.m_evo)
R_CO10_CI = 0.18*77.83 # r=0.18 times cubic frequency scaling
xRJ_GHz = 0.0176086

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

bandpass_PA6_090_edges = [bandpass_PA6_090[1,1].-diff(bandpass_PA6_090[1:2,1])/2;bandpass_PA6_090[1:end-1,1]+diff(bandpass_PA6_090[:,1])/2;bandpass_PA6_090[end,1].+diff(bandpass_PA6_090[end-1:end,1])/2]
bandpass_PA6_150_edges = [bandpass_PA6_150[1,1].-diff(bandpass_PA6_150[1:2,1])/2;bandpass_PA6_150[1:end-1,1]+diff(bandpass_PA6_150[:,1])/2;bandpass_PA6_150[end,1].+diff(bandpass_PA6_150[end-1:end,1])/2]
bandpass_PA4_220_edges = [bandpass_PA4_220[1,1].-diff(bandpass_PA4_220[1:2,1])/2;bandpass_PA4_220[1:end-1,1]+diff(bandpass_PA4_220[:,1])/2;bandpass_PA4_220[end,1].+diff(bandpass_PA4_220[end-1:end,1])/2]

bandpass_PA6_090_norm = bandpass_PA6_090[:,2]/sum(bandpass_PA6_090[:,2])
bandpass_PA6_150_norm = bandpass_PA6_150[:,2]/sum(bandpass_PA6_150[:,2])
bandpass_PA4_220_norm = bandpass_PA4_220[:,2]/sum(bandpass_PA4_220[:,2])

bandpass_PA6_090_diff = diff(bandpass_PA6_090_edges)
bandpass_PA6_150_diff = diff(bandpass_PA6_150_edges)
bandpass_PA4_220_diff = diff(bandpass_PA4_220_edges)

bandpass_PA6_090_len = length(bandpass_PA6_090_norm)
bandpass_PA6_150_len = length(bandpass_PA6_150_norm)
bandpass_PA4_220_len = length(bandpass_PA4_220_norm)

function write_chunk(
                     output_dir, chunk_index, model, cosmo,
                     pos, mass, freqs)
    sources = generate_sources(model, cosmo, pos, mass);
    fluxes_cen = Array{Float32,1}(undef,sources.N_cen)
    fluxes_sat = Array{Float32,1}(undef,sources.N_sat)
    LIR_cen = Array{Float32,1}(undef,sources.N_cen)
    LIR_sat = Array{Float32,1}(undef,sources.N_sat)
    LcoJ_cen = Array{Float32,2}(undef,sources.N_cen,7)
    LcoJ_sat = Array{Float32,2}(undef,sources.N_sat,7)
    Lnu_to_LIR(Td) = 9.521f8*(53.865 + 3.086*Td - 0.213*Td^2 + 0.00749*Td^3);
    # calculate observing frequencies and quasi-temperatures for each halo
    nuJ_cen = Array{Float32,2}(undef,(sources.N_cen,7))
    nuJ_sat = Array{Float32,2}(undef,(sources.N_sat,7))
    quasiTcoJ_cen = Array{Float32,2}(undef,(sources.N_cen,7))
    quasiTcoJ_sat = Array{Float32,2}(undef,(sources.N_sat,7))
    fluxCO_090_cen = Array{Float32,2}(undef,(sources.N_cen,7))
    fluxCO_090_sat = Array{Float32,2}(undef,(sources.N_sat,7))
    fill!(fluxCO_090_cen,zero(Float32))
    fill!(fluxCO_090_sat,zero(Float32))
    fluxCO_150_cen = Array{Float32,2}(undef,(sources.N_cen,7))
    fluxCO_150_sat = Array{Float32,2}(undef,(sources.N_sat,7))
    fill!(fluxCO_150_cen,zero(Float32))
    fill!(fluxCO_150_sat,zero(Float32))
    fluxCO_220_cen = Array{Float32,2}(undef,(sources.N_cen,7))
    fluxCO_220_sat = Array{Float32,2}(undef,(sources.N_sat,7))
    fill!(fluxCO_220_cen,zero(Float32))
    fill!(fluxCO_220_sat,zero(Float32))
    nuCI_cen = Array{Float32,1}(undef,sources.N_cen)
    nuCI_sat = Array{Float32,1}(undef,sources.N_sat)
    LCI_cen = Array{Float32,1}(undef,sources.N_cen)
    LCI_sat = Array{Float32,1}(undef,sources.N_sat)
    quasiTCI_cen = Array{Float32,1}(undef,sources.N_cen)
    quasiTCI_sat = Array{Float32,1}(undef,sources.N_sat)
    fluxCI_090_cen = Array{Float32,1}(undef,sources.N_cen)
    fill!(fluxCI_090_cen,zero(Float32))
    fluxCI_150_cen = Array{Float32,1}(undef,sources.N_cen)
    fill!(fluxCI_150_cen,zero(Float32))
    fluxCI_220_cen = Array{Float32,1}(undef,sources.N_cen)
    fill!(fluxCI_220_cen,zero(Float32))
    fluxCI_090_sat = Array{Float32,1}(undef,sources.N_sat)
    fill!(fluxCI_090_sat,zero(Float32))
    fluxCI_150_sat = Array{Float32,1}(undef,sources.N_sat)
    fill!(fluxCI_150_sat,zero(Float32))
    fluxCI_220_sat = Array{Float32,1}(undef,sources.N_sat)
    fill!(fluxCI_220_sat,zero(Float32))
    Threads.@threads :static for i in 1:sources.N_cen
        # calculate dust temperatures and bolometric LIR for each source
        Td = model.shang_Td * (1.f0 .+sources.redshift_cen[i])^model.shang_alpha;
        LIR_cen[i] = Lnu_to_LIR(Td)*sources.lum_cen[i]*XGPaint.nu2theta(2.10833f12,sources.redshift_cen[i],model)
        rnd = Float32(randn())
        r21 = 0.75f0+0.11f0*(rnd/2.0f0+373.f0/275.f0-sqrt(rnd^2.0f0/4.0f0-252.0f0/275.0f0*rnd+139129.0f0/75625.0f0))
        # generate CO SLED for each source with log-normal LCO(1-0) scatter
        LcoJ_cen[i,1] = LIR_cen[i]^(1/1.37)*10^(1.74/1.37)*exp.((randn().-0.5*2.302585*0.3)*2.302585*0.3)*4.9e-5
        for J in 2:7
            LcoJ_cen[i,J] = LcoJ_cen[i,1]/(1.0f0+exp((log.(1.0f0/r21-1.0f0)-1.0f0)+Float32(J-1)))*J^3
        end
        for J in 1:7
            nuJ_cen[i,J] = 115.27f0*J/(1.0f0+sources.redshift_cen[i])
            quasiTcoJ_cen[i,J] = LcoJ_cen[i,J]*(1.315415f0/4.0f0/pi/nuJ_cen[i,J]^2.0f0/sources.dist_cen[i]^2.0f0/(1.0f0+sources.redshift_cen[i])^2.0f0) # divide by dnu in GHz
            xRJ = xRJ_GHz*nuJ_cen[i,J]
            quasiTcoJ_cen[i,J]/= xRJ^2*exp(xRJ)/(exp(xRJ)-1)^2
        end
        nuCI_cen[i] = 492.16f0/(1.0f0+sources.redshift_cen[i])
        LCI_cen[i] = LcoJ_cen[i,1]*R_CO10_CI*exp((randn()-0.5*2.302585*0.2)*2.302585*0.2)
        quasiTCI_cen[i] = LCI_cen[i]*(1.315415f0/4.0f0/pi/nuCI_cen[i]^2.0f0/sources.dist_cen[i]^2.0f0/(1.0f0+sources.redshift_cen[i])^2.0f0) # divide by dnu in GHz
        xRJ = xRJ_GHz*nuCI_cen[i]
        quasiTCI_cen[i]/= xRJ^2*exp(xRJ)/(exp(xRJ)-1)^2
    end
    Threads.@threads :static for i in 1:sources.N_sat
        Td = model.shang_Td * (1.f0 .+sources.redshift_sat[i])^model.shang_alpha;
        LIR_sat[i] = Lnu_to_LIR(Td)*sources.lum_sat[i]*XGPaint.nu2theta(2.10833f12,sources.redshift_sat[i],model)
        rnd = Float32(randn())
        r21 = 0.75f0+0.11f0*(rnd/2.0f0+373.f0/275.f0-sqrt(rnd^2.0f0/4.0f0-252.0f0/275.0f0*rnd+139129.0f0/75625.0f0))
        LcoJ_sat[i,1] = LIR_sat[i]^(1/1.37)*10^(1.74/1.37)*exp((randn()-0.5*2.302585*0.3)*2.302585*0.3)*4.9e-5
        for J in 2:7
            LcoJ_sat[i,J] = LcoJ_sat[i,1]/(1.0f0+exp((log.(1.0f0/r21-1.0f0)-1.0f0)+Float32(J-1)))*J^3
        end
        for J in 1:7
            nuJ_sat[i,J] = 115.27f0*J/(1.0f0+sources.redshift_sat[i])
            quasiTcoJ_sat[i,J] = LcoJ_sat[i,J]*(1.315415f0/4.0f0/pi/nuJ_sat[i,J]^2.0f0/sources.dist_sat[i]^2.0f0/(1.0f0+sources.redshift_sat[i])^2.0f0) # divide by dnu in GHz
            xRJ = xRJ_GHz*nuJ_sat[i,J]
            quasiTcoJ_sat[i,J]/= xRJ^2*exp(xRJ)/(exp(xRJ)-1)^2
        end
        nuCI_sat[i] = 492.16f0/(1.0f0+sources.redshift_sat[i])
        LCI_sat[i] = LcoJ_sat[i,1]*R_CO10_CI*exp((randn()-0.5*2.302585*0.2)*2.302585*0.2)
        quasiTCI_sat[i] = LCI_sat[i]*(1.315415f0/4.0f0/pi/nuCI_sat[i]^2.0f0/sources.dist_sat[i]^2.0f0/(1.0f0+sources.redshift_sat[i])^2.0f0) # divide by dnu in GHz
        xRJ = xRJ_GHz*nuCI_sat[i]
        quasiTCI_sat[i]/= xRJ^2*exp(xRJ)/(exp(xRJ)-1)^2
    end
    println("writing info for ",sources.N_cen," centrals")
    h5open(joinpath(output_dir, "sources/cen_chunk$(chunk_index).h5"), "w") do file
        write(file, "redshift", sources.redshift_cen)
        write(file, "theta", sources.theta_cen)
        write(file, "phi", sources.phi_cen)
        write(file, "LIR", LIR_cen)
        write(file, "LCO", LcoJ_cen)
        write(file, "LCI", LCI_cen)
    end
    println("writing info for ",sources.N_sat," satellites")
    h5open(joinpath(output_dir, "sources/sat_chunk$(chunk_index).h5"), "w") do file
        write(file, "redshift", sources.redshift_sat)
        write(file, "theta", sources.theta_sat)
        write(file, "phi", sources.phi_sat)
        write(file, "LIR", LIR_sat)
        write(file, "LCO", LcoJ_sat)
        write(file, "LCI", LCI_sat)
    end
    m = HealpixMap{Float64,RingOrder}(model.nside)
    mCO090 = HealpixMap{Float64,RingOrder}(model.nside)
    mCO150 = HealpixMap{Float64,RingOrder}(model.nside)
    mCO220 = HealpixMap{Float64,RingOrder}(model.nside)
    mCO090_pixsr = nside2pixarea(mCO090.resolution.nside)
    mCO150_pixsr = nside2pixarea(mCO150.resolution.nside)
    mCO220_pixsr = nside2pixarea(mCO220.resolution.nside)
    mCI090 = HealpixMap{Float64,RingOrder}(model.nside)
    mCI150 = HealpixMap{Float64,RingOrder}(model.nside)
    mCI220 = HealpixMap{Float64,RingOrder}(model.nside)
    println("painting CO/CI from ",sources.N_cen," centrals")
    Threads.@threads :static for i in 1:sources.N_cen
        for J in 1:7
            j = searchsortedlast(bandpass_PA6_090_edges,nuJ_cen[i,J])
            if (j>0) && (j<=bandpass_PA6_090_len)
                fluxCO_090_cen[i,J] = quasiTcoJ_cen[i,J]*bandpass_PA6_090_norm[j]/bandpass_PA6_090_diff[j]/mCO090_pixsr;
                mCO090[sources.hp_ind_cen[i]]+= fluxCO_090_cen[i,J];
            end
            j = searchsortedlast(bandpass_PA6_150_edges,nuJ_cen[i,J])
            if (j>0) && (j<=bandpass_PA6_150_len)
                fluxCO_150_cen[i,J] = quasiTcoJ_cen[i,J]*bandpass_PA6_150_norm[j]/bandpass_PA6_150_diff[j]/mCO150_pixsr;
                mCO150[sources.hp_ind_cen[i]]+= fluxCO_150_cen[i,J];
            end
            j = searchsortedlast(bandpass_PA4_220_edges,nuJ_cen[i,J])
            if (j>0) && (j<=bandpass_PA4_220_len)
                fluxCO_220_cen[i,J] = quasiTcoJ_cen[i,J]*bandpass_PA4_220_norm[j]/bandpass_PA4_220_diff[j]/mCO220_pixsr;
                mCO220[sources.hp_ind_cen[i]]+= fluxCO_220_cen[i,J];
            end
        end
        j = searchsortedlast(bandpass_PA6_090_edges,nuCI_cen[i])
        if (j>0) && (j<=bandpass_PA6_090_len)
            fluxCI_090_cen[i] = quasiTCI_cen[i]*bandpass_PA6_090_norm[j]/bandpass_PA6_090_diff[j]/mCO090_pixsr;
            mCI090[sources.hp_ind_cen[i]]+= fluxCI_090_cen[i];
        end
        j = searchsortedlast(bandpass_PA6_150_edges,nuCI_cen[i])
        if (j>0) && (j<=bandpass_PA6_150_len)
            fluxCI_150_cen[i] = quasiTCI_cen[i]*bandpass_PA6_150_norm[j]/bandpass_PA6_150_diff[j]/mCO150_pixsr;
            mCI150[sources.hp_ind_cen[i]]+= fluxCI_150_cen[i];
        end
        j = searchsortedlast(bandpass_PA4_220_edges,nuCI_cen[i])
        if (j>0) && (j<=bandpass_PA4_220_len)
            fluxCI_220_cen[i] = quasiTCI_cen[i]*bandpass_PA4_220_norm[j]/bandpass_PA4_220_diff[j]/mCO220_pixsr;
            mCI220[sources.hp_ind_cen[i]]+= fluxCI_220_cen[i];
        end
    end
    println("painting CO/CI from ",sources.N_sat," satellites")
    Threads.@threads :static for i in 1:sources.N_sat
        for J in 1:7
            j = searchsortedlast(bandpass_PA6_090_edges,nuJ_sat[i,J])
            if (j>0) && (j<=bandpass_PA6_090_len)
                fluxCO_090_sat[i,J] = quasiTcoJ_sat[i,J]*bandpass_PA6_090_norm[j]/bandpass_PA6_090_diff[j]/mCO090_pixsr;
                mCO090[sources.hp_ind_sat[i]]+= fluxCO_090_sat[i,J];
            end
            j = searchsortedlast(bandpass_PA6_150_edges,nuJ_sat[i,J])
            if (j>0) && (j<=bandpass_PA6_150_len)
                fluxCO_150_sat[i,J] = quasiTcoJ_sat[i,J]*bandpass_PA6_150_norm[j]/bandpass_PA6_150_diff[j]/mCO150_pixsr;
                mCO150[sources.hp_ind_sat[i]]+= fluxCO_150_sat[i,J];
            end
            j = searchsortedlast(bandpass_PA4_220_edges,nuJ_sat[i,J])
            if (j>0) && (j<=bandpass_PA4_220_len)
                fluxCO_220_sat[i,J] = quasiTcoJ_sat[i,J]*bandpass_PA4_220_norm[j]/bandpass_PA4_220_diff[j]/mCO220_pixsr;
                mCO220[sources.hp_ind_sat[i]]+= fluxCO_220_sat[i,J];
            end
        end
        j = searchsortedlast(bandpass_PA6_090_edges,nuCI_sat[i])
        if (j>0) && (j<=bandpass_PA6_090_len)
            fluxCI_090_sat[i] = quasiTCI_sat[i]*bandpass_PA6_090_norm[j]/bandpass_PA6_090_diff[j]/mCO090_pixsr;
            mCI090[sources.hp_ind_sat[i]]+= fluxCI_090_sat[i];
        end
        j = searchsortedlast(bandpass_PA6_150_edges,nuCI_sat[i])
        if (j>0) && (j<=bandpass_PA6_150_len)
            fluxCI_150_sat[i] = quasiTCI_sat[i]*bandpass_PA6_150_norm[j]/bandpass_PA6_150_diff[j]/mCO150_pixsr;
            mCI150[sources.hp_ind_sat[i]]+= fluxCI_150_sat[i];
        end
        j = searchsortedlast(bandpass_PA4_220_edges,nuCI_sat[i])
        if (j>0) && (j<=bandpass_PA4_220_len)
            fluxCI_220_sat[i] = quasiTCI_sat[i]*bandpass_PA4_220_norm[j]/bandpass_PA4_220_diff[j]/mCO220_pixsr;
            mCI220[sources.hp_ind_sat[i]]+= fluxCI_220_sat[i];
        end
    end
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
            write(file, "flux", fluxCO_090_cen[:,J]; compress=true)
        end
        jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCO$(J)_150.jld2"), "w") do file
            write(file, "flux", fluxCO_150_cen[:,J]; compress=true)
        end
        jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCO$(J)_220.jld2"), "w") do file
            write(file, "flux", fluxCO_220_cen[:,J]; compress=true)
        end
        jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCO$(J)_090.jld2"), "w") do file
            write(file, "flux", fluxCO_090_sat[:,J]; compress=true)
        end
        jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCO$(J)_150.jld2"), "w") do file
            write(file, "flux", fluxCO_150_sat[:,J]; compress=true)
        end
        jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCO$(J)_220.jld2"), "w") do file
            write(file, "flux", fluxCO_220_sat[:,J]; compress=true)
        end
    end
    println("writing CI fluxes ...")
    jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCI_090.jld2"), "w") do file
        write(file, "flux", fluxCI_090_cen; compress=true)
    end
    jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCI_150.jld2"), "w") do file
        write(file, "flux", fluxCI_150_cen; compress=true)
    end
    jldopen(joinpath(output_dir, "sources/cen_chunk$(chunk_index)_fluxCI_220.jld2"), "w") do file
        write(file, "flux", fluxCI_220_cen; compress=true)
    end
    jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCI_090.jld2"), "w") do file
        write(file, "flux", fluxCI_090_sat; compress=true)
    end
    jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCI_150.jld2"), "w") do file
        write(file, "flux", fluxCI_150_sat; compress=true)
    end
    jldopen(joinpath(output_dir, "sources/sat_chunk$(chunk_index)_fluxCI_220.jld2"), "w") do file
        write(file, "flux", fluxCI_220_sat; compress=true)
    end
    for freq in freqs

        XGPaint.paint!(m, parse(Float32, freq) * 1.0f9, model, sources,
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

## the rest of this is more or less the same as cib_planck2013_chunked.jl
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
        write_chunk(output_dir, chunk_index, model, cosmo,
            pos, mass, freqs)
    end
    # save maps
    println("writing CO+CI maps ...")
    m = HealpixMap{Float64,RingOrder}(model.nside)
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
scratch_dir = "/home/dongwooc/scratchspace/cib_co_sources_scarfy"
println("SCRATCH: ", scratch_dir)
mkpath(scratch_dir)
mkpath(joinpath(scratch_dir, "sources"))
run_all_chunks(scratch_dir, halo_pos, halo_mass, freqs)
