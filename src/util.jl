using HDF5
using Healpix
using Random
using Random: MersenneTwister
using LazyArtifacts

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


function getrange(n)
    tid = Threads.threadid()
    nt = Threads.nthreads()
    d , r = divrem(n, nt)
    from = (tid - 1) * d + min(r, tid - 1) + 1
    to = from + d - 1 + (tid ≤ r ? 1 : 0)
    from:to
end


function threaded_rand!(random_number_generators, arr::Array{T,1};
      chunksize=4096) where T

   num = size(arr,1)
   Threads.@threads :static for (i1, i2) in chunk(num, chunksize)
      @views rand!(random_number_generators[Threads.threadid()], arr[i1:i2])
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


@doc raw"""
    queryDiscRing!(result, ringinfo, resol::Resolution, theta, phi, radius; fact=0)

In-place calculation of a list of the indices of those pixels whose centers are closer
than `radius` to direction `(theta, phi)`. The three angles `radius`,
`theta`, and `phi` must be expressed in radians.

If `fact` is nonzero, it must be a positive integer; it requires to
carry the computation at a resolution `fact * nside`.
"""
function queryDiscRing!(
    result::AbstractArray{Int,1},
    ringinfo,
    resol::Resolution,
    theta,
    phi,
    radius;
    fact=0,
)
    inclusive = (fact != 0)
    empty!(result)

    fct = 1

    if inclusive
        @assert ((1 << ORDER_MAX) / resol.nside) >= fact
        fct = fact
    end

    b2 = Resolution(fct * resol.nside)
    (rsmall, rbig) = if fct > 1
        (radius + max_pixrad(b2), radius + max_pixrad(resol))
    else
        value = inclusive ? (radius + max_pixrad(resol)) : radius
        (value, value)
    end
    
    (rsmall >= π) && (return 1:resol.npix)

    rbig = min(pi, rbig)
    (cosrsmall, cosrbig) = (cos(rsmall), cos(rbig))

    z0 = cos(theta)
    xa = 1 / sqrt((1 - z0) * (1 + z0))

    cpix = zphi2pixRing(resol, z0, phi)

    rlat1 = theta - rsmall
    zmax = cos(rlat1)
    irmin = Healpix.ringAbove(resol, zmax) + 1

    if (rlat1 <= 0) && (irmin > 1)
        getringinfo!(resol, irmin - 1, ringinfo)
        append!(
            result,
            ringinfo.firstPixIdx:(ringinfo.firstPixIdx + ringinfo.numOfPixels - 1),
        )
    end

    if (fct > 1) && (rlat1 > 0)
        irmin = max(1, irmin - 1)
    end

    rlat2 = theta + rsmall
    zmin = cos(rlat2)
    irmax = Healpix.ringabove(resol, zmin)

    if (fct > 1) && (rlat2 < π)
        irmax = min(resol.nsideTimesFour - 1, irmax + 1)
    end

    for iz in irmin:irmax
        z = ring2z(resol, iz)
        x = (cosrbig - z * z0) * xa
        ysq = 1 - z^2 - x^2
        dphi = if ysq < 0
            (fct == 1) ? 0 : π - 1e-15
        else
            atan(sqrt(ysq), x)
        end

        if dphi > 0
            getringinfo!(resol, iz, ringinfo)
            ipix1 = ringinfo.firstPixIdx
            nr = ringinfo.numOfPixels
            shift = ringinfo.shifted ? 0.5 : 0.0

            ipix2 = ipix1 + nr - 1

            ip_lo = floor(Int, nr / 2π * (phi - dphi) - shift) + 1
            ip_hi = floor(Int, nr / 2π * (phi + dphi) - shift)

            if fct > 1
                while ((ip_lo <= ip_hi) &&
                       checkPixelRing(resol, b2, ip_lo, nr, ipix1, fct, z0, phi, cosrsmall, cpix))
                    ip_lo += 1
                end
                while ((ip_hi > ip_lo) &&
                       checkPixelRing(resol, b2, ip_hi, nr, ipix1, fct, z0, phi, cosrsmall, cpix))
                    ip_hi -= 1
                end
            end

            if ip_lo <= ip_hi
                if ip_hi >= nr
                    ip_lo -= nr
                    ip_hi -= nr
                end

                if ip_lo < 0
                    append!(result, ipix1:(ipix1 + ip_hi))
                    append!(result, (ipix1 + ip_lo + nr):ipix2)
                else
                    append!(result, (ipix1 + ip_lo):(ipix1 + ip_hi))
                end
            end
        end

        if (rlat2 >= π) && (irmax + 1 < resol.nsideTimesFour)
            getringinfo!(resol, irmax + 1, ringinfo)
            append!(result, ringinfo.firstPixIdx:(ringinfo.firstPixIdx + ringinfo.numOfPixels - 1))
        end
    end

    result
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

