using HDF5
using Healpix
using Random
using Random: MersenneTwister
using LazyArtifacts
using ChunkSplitters: chunks

"""
Utility function to read an HDF5 table with x, y, z, M_h as the four rows.
The hdf5 record is "halos".
"""
function read_halo_catalog_hdf5(filename)
    hdata = h5open(filename, "r") do file
        read(file, "halos")
    end
    pos = hdata[1:3,:]
    halo_mass = hdata[4,:]
    return pos, halo_mass
end

"""
Reads a collection of example halos out of the WebSky halo catalogs, and returns 
RA (rad), DEC (rad), redshift, and halo mass (M200c).
"""
function load_example_halos()
    rootpath = artifact"tsz_example"
    fid = h5open(joinpath(rootpath, "tsz_example", "little_box_m200c_v2.h5"), "r")
    ra, dec = deg2rad.(fid["ra"]), deg2rad.(fid["dec"])
    redshift, halo_mass = collect(fid["redshift"]), collect(fid["halo_mass"])
    close(fid)

    return ra, dec, redshift, halo_mass
end

function load_example_tsz_map()
    rootpath = artifact"tsz_test_map"
    ymap = jldopen(joinpath(rootpath, "test_tsz_map.jld2"), "r")["data"]
    return ymap
end


"""
Returns a stored artifact: a precomputed SZpack table
"""
function rsz_szpack_table_filename()
    return joinpath(artifact"rsz_table", "szpack_interp_T75_upd.dat")
end


function load_precomputed_battaglia_data()
    rootpath = artifact"tsz_example"
    model_file = joinpath(rootpath, "tsz_example", "battaglia_interpolation.jld2")
    model = load(model_file)
    prof_logθs, prof_redshift, prof_logMs, prof_y = model["prof_logθs"], 
        model["prof_redshift"], model["prof_logMs"], model["prof_y"]
    return prof_y, prof_logθs, prof_redshift, prof_logMs
end


"""
Reads in a standard Battaglia 2016 model from disk for tSZ.
"""
function load_precomputed_battaglia()
    prof_y, prof_logθs, prof_redshift, prof_logMs = load_precomputed_battaglia_data()
    interp_model = scale(Interpolations.interpolate(log.(prof_y), BSpline(Cubic(Line(OnGrid())))), 
        prof_logθs, prof_redshift, prof_logMs);
    p = Battaglia16ThermalSZProfile(Omega_c=0.2589, Omega_b=0.0486, h=0.6774)
    return LogInterpolatorProfile(p, interp_model)
end


function load_precomputed_battaglia_tau_data()
    rootpath = artifact"precomputed_battaglia_tau"
    model_file = joinpath(rootpath, "battaglia_tau_profile_interp.jld2")
    model = load(model_file)
    prof_logθs, prof_redshift, prof_logMs, prof_y = model["prof_logθs"], 
        model["prof_redshift"], model["prof_logMs"], model["prof_y"]
    return prof_y, prof_logθs, prof_redshift, prof_logMs
end


function load_precomputed_battaglia_tau()
    prof_y, prof_logθs, prof_redshift, prof_logMs = load_precomputed_battaglia_tau_data()
    interp_model = scale(Interpolations.interpolate(log.(prof_y), BSpline(Cubic(Line(OnGrid())))), 
        prof_logθs, prof_redshift, prof_logMs);
    p = BattagliaTauProfile(Omega_c=0.2589, Omega_b=0.0486, h=0.6774)
    return LogInterpolatorProfile(p, interp_model)
end

"""
Generates a list of tuples which describe starting and ending chunk indices.
Useful for parallelizing an array operation.
"""
function chunk(arr_len, chunksize::Integer)
    return [(i, min(i + chunksize-1, arr_len))
        for i in range(1, arr_len, step=chunksize)]
end


function threaded_rand!(random_number_generators, arr::Array{T,1}) where T

   # Use ChunkSplitters for better threading
   Threads.@threads for chunk in chunks(eachindex(arr); n=Threads.nthreads())
      tid = Threads.threadid()
      @views rand!(random_number_generators[tid], arr[chunk])
   end
end

# if no RNG is supplied, make some
function threaded_rand!(arr::Array{T,1}; chunksize=4096) where T
    random_number_generators = [Random.default_rng(i) for i in 1: Threads.nthreads()]
    threaded_rand!(random_number_generators, arr; chunksize=chunksize)
end

"""
Generate an array where the value at index i corresponds to the index of the
first source of halo i. Takes an array where the value at index i corresponds
to the number of subhalos that halo i has.
"""
function generate_subhalo_offsets(num_subhalos)
    result = cumsum(num_subhalos)
    prepend!(result, 0)
    return result
end

"""
Utility function which prepends some zeros to an array. It makes a copy instead
of modifying the input.
"""
function ellpad(arr::Array{T,N}; nzeros=1) where {T,N}
    result = arr[:]
    pushfirst!(result, zeros(T, nzeros)...)
    return result
end

function sort_halo_catalog(ra, dec, redshift, halo_mass)
    perm = sortperm(dec)  # sortperm(dec, alg=ThreadsX.MergeSort)
    ra = ra[perm]
    dec = dec[perm]
    redshift = redshift[perm]
    halo_mass = halo_mass[perm]
    return ra, dec, redshift, halo_mass
end


function catalog2map!(m::HealpixMap{T,RingOrder}, flux, theta, phi) where T
    res = m.resolution
    pixel_array = m.pixels
    N_halo = length(flux)
    per_pixel_steradian = 1 / nside2pixarea(res.nside)  # divide by healpix pixel size

    # try to prevent thread issues by sorting by theta
    perm = sortperm(theta, rev=true, alg=ThreadsX.MergeSort)
    Threads.@threads :static for i_perm in 1:N_halo
        i_halo = perm[i_perm]
        hp_ind = Healpix.ang2pixRing(res, theta[i_halo], phi[i_halo])
        pixel_array[hp_ind] += flux[i_halo] * per_pixel_steradian
    end
end


function catalog2map!(m::Enmap{T}, flux, theta, phi, pixsizes; erase_first=true) where T
    N_halo = length(flux)
    pixel_array = parent(m)
    
    if erase_first
        fill!(pixel_array, zero(T))
    end

    α = T.(phi)
    δ = T(π/2) .- T.(theta) 

    ipix, jpix = sky2pix(size(m), m.wcs, α, δ)

    len1, len2 = size(m)

    # try to prevent thread issues by sorting by theta
    perm = sortperm(theta, rev=true, alg=ThreadsX.MergeSort)
    Threads.@threads :static for i_perm in 1:N_halo
        i_halo = perm[i_perm]
        i = round(Int, ipix[i_halo])
        i = mod(i - 1, size(m, 1)) + 1
        j = round(Int, jpix[i_halo])
        if (1 ≤ i ≤ len1) && (1 ≤ j ≤ len2)
            pixel_array[i, j] += flux[i_halo] / pixsizes[i, j]
        end
    end
    
    return m
end



"""
Compute the real-space beam from a harmonic-space beam.
"""
function bl2beam(bl::AbstractArray{T,1}, θ::AbstractArray) where T
    lmax = length(bl) - 1
    nx = length(θ)
    x = cos.(θ)
    p0 = ones(T, nx)
    p1 = x

    beam = bl[0+1] * p0 + bl[1+1] * p1 * 3

    for ℓ in 2:lmax
        p2 = @. x * p1 * (2ℓ - 1) / ℓ - p0 * (ℓ - 1) / ℓ
        p0 = p1
        p1 = p2
        beam += bl[ℓ+1] * p2 * (2ℓ + 1)
    end

    beam /= 4 * T(π)
    return beam
end


"""
    vecmap([T=Float64], nside::Int)

Generate a map of 3-tuples, showing the (x,y,z) values of the points on the unit 2-sphere
for a Healpix map.
"""
function vectorhealpixmap(::Type{T}, nside::Int) where T
    npix = nside2npix(nside)
    res = Resolution(nside)
    arr = Vector{Tuple{T, T, T}}(undef, npix)
    posmap = HealpixMap{Tuple{T, T, T}, RingOrder}(arr)
    Threads.@threads :static for i in 1:npix
        posmap.pixels[i] = pix2vecRing(res, i)
    end
    return posmap.pixels
end
vectorhealpixmap(nside::Int) = vectorhealpixmap(Float64, nside)


export read_halo_catalog_hdf5, chunk, generate_subhalo_offsets

