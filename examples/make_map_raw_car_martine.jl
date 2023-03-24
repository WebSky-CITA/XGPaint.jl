# modified from Zack Li

using Pixell, WCS, XGPaint
using Cosmology
using Interpolations
import XGPaint: AbstractProfile
using HDF5
import JSON
using JLD2
using ThreadsX
using FileIO
using Healpix

print("Threads: ", Threads.nthreads(), "\n")
modeltype::String = "break"
healpix = true


if healpix
    nside = 1024
    m = HealpixMap{Float64,RingOrder}(nside)
else
    shape, wcs = fullsky_geometry(0.5 * Pixell.arcminute)
    # for a box:
    # box = [30   -30;           # RA
    #        -25     25] * Pixell.degree  # DEC
    # shape, wcs = geometry(CarClenshawCurtis{Float64}, box, 0.5 * Pixell.arcminute)
    m = Enmap(zeros(shape), wcs)
    # precomputed sky angles
    α_map, δ_map = posmap(shape, wcs)
    psa = (sin_α=sin.(α_map), cos_α=cos.(α_map), sin_δ=sin.(δ_map), cos_δ=cos.(δ_map))
end

##
fid = h5open("/mnt/raid-cita/mlokken/buzzard/catalogs/halos/buzzard_halos.hdf5", "r")
# fid = h5open("/mnt/scratch-lustre/mlokken/pkpatch/halos_hdf5.h5", "r")
ra, dec = deg2rad.(fid["ra"]), deg2rad.(fid["dec"])
redshift = collect(fid["z"])
halo_mass = collect(fid["m200c"])
# choose the cutoff for the pressure profiles here
cutoff = 4 # the standard
apply_beam = true
##
perm = sortperm(dec, alg=ThreadsX.MergeSort)
ra = ra[perm]
dec = dec[perm]
redshift = redshift[perm]
halo_mass = halo_mass[perm]



##
print("Precomputing the model profile grid.\n")

# set up a profile to paint
if modeltype=="break"
    p = XGPaint.BreakModel(Omega_c=0.24, Omega_b=0.046, h=0.7, alpha_break=1.5) # Buzzard cosmology
elseif modeltype=="battaglia"
    p = XGPaint.BattagliaProfile(Omega_c=0.24, Omega_b=0.046, h=0.7) # Buzzard cosmology
end

# beam stuff
N_logθ = 512
rft = RadialFourierTransform(n=N_logθ, pad=256)
beamstr = ""

if apply_beam
    ells = rft.l
    fwhm = deg2rad(1.6/60.)  # ACT beam
    σ² = fwhm^2 / (8log(2))
    lbeam = @. exp(-ells * (ells+1) * σ² / 2)
    beamstr = "_1p6arcmin"
end
model_file::String = "cached_$(modeltype)_$(beamstr).jld2"
if isfile(model_file)
    print("Found cached Battaglia profile model. Loading from disk.\n")
    model = load(model_file)
    prof_logθs, prof_redshift, prof_logMs, prof_y = model["prof_logθs"], 
        model["prof_redshift"], model["prof_logMs"], model["prof_y"]
else
    print("Didn't find a cached profile model. Computing and saving.\n")
    logθ_min, logθ_max = log(minimum(rft.r)), log(maximum(rft.r))
    @time prof_logθs, prof_redshift, prof_logMs, prof_y = profile_grid(p; 
        N_logθ=length(rft.r), logθ_min=logθ_min, logθ_max=logθ_max)
    
    if apply_beam
        # now apply the beam
        @time XGPaint.transform_profile_grid!(prof_y, rft, lbeam)
        prof_y = prof_y[begin+rft.pad:end-rft.pad, :, :]
        prof_logθs = prof_logθs[begin+rft.pad:end-rft.pad]
        XGPaint.cleanup_negatives!(prof_y)
    end
    jldsave(model_file; prof_logθs, prof_redshift, prof_logMs, prof_y)
    
end

if ~healpix # do this for CAR
    itp = Interpolations.interpolate(log.(prof_y), BSpline(Cubic(Line(OnGrid()))))
    sitp = scale(itp, prof_logθs, prof_redshift, prof_logMs);
end


##
function paint_map!(m::Enmap, p::XGPaint.AbstractProfile, psa, sitp, masses, 
                    redshifts, αs, δs, irange; mult=4)
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax = XGPaint.θmax(p, mh * XGPaint.M_sun, z, mult=mult)
        profile_paint!(m, α₀, δ₀, psa, sitp, z, mh, θmax)
    end
end

function paint_map!(m::HealpixMap, p::XGPaint.AbstractProfile{T}, masses, 
    redshifts, αs, δs, irange; mult=4) where {T}
    for i in irange
        α₀ = αs[i]
        δ₀ = δs[i]
        mh = masses[i]
        z = redshifts[i]
        θmax = XGPaint.θmax(p, mh * XGPaint.M_sun, z, mult=mult)
        θmin = 0.0 # don't know what this should be # error
        beam_interp = XGPaint.realspacegaussbeam(T,fwhm)[0] #error
        w = XGPaint.HealpixPaintingWorkspace(nside, θmin, θmax, beam_interp)
        profile_paint!(m, α₀, δ₀, w, z, mh, θmax)    
    end
end

function chunked_paint!(m::Enmap, p::XGPaint.AbstractProfile, psa, sitp, masses, 
                        redshifts, αs, δs; mult=4)
    m .= 0.0
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_map!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2, mult=mult)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_map!(m, p, psa, sitp, masses, redshifts, αs, δs, i1:i2, mult=mult)
    end
end

function chunked_paint!(m::HealpixMap, p::XGPaint.AbstractProfile, masses, 
                        redshifts, αs, δs; mult=4)
    m .= 0.0
    
    N_sources = length(masses)
    chunksize = ceil(Int, N_sources / (2Threads.nthreads()))
    chunks = chunk(N_sources, chunksize);
    
    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i
        i1, i2 = chunks[chunk_i]
        paint_map!(m, p, masses, redshifts, αs, δs, i1:i2, mult=mult)
    end

    Threads.@threads for i in 1:Threads.nthreads()
        chunk_i = 2i - 1
        i1, i2 = chunks[chunk_i]
        paint_map!(m, p,  masses, redshifts, αs, δs, i1:i2, mult=mult)
    end
end
##
# cut = eachindex(halo_mass)[begin:end-5]

print(maximum(ra), minimum(ra), maximum(dec), minimum(dec))
print("Painting map.\n")
if healpix
    @time chunked_paint!(m, p, halo_mass, redshift, ra, dec, mult=cutoff) # for a Healpix
else
    @time chunked_paint!(m, p, psa, sitp, halo_mass, redshift, ra, dec, mult=cutoff) # for an Enmap
end

#
write_map(
    "/mnt/raid-cita/mlokken/buzzard/ymaps/ymap_buzzard_$(modeltype)_bbps_car$(beamstr)_cutoff$(cutoff).fits",
    m)

using PyPlot
plt.clf()
plt.figure()
plt.imshow(log10.(m.data'))
plt.axis("off")
plt.savefig("test_Healpix.png", bbox_inches="tight",pad_inches = 0)
plt.gcf()



# ##
# m .= 0
# profile_paint!(m, 0., 0., psa, sitp, 0.1, 1e14, π/90)

# ##
# using Plots
# Plots.plot(log10.(m))
