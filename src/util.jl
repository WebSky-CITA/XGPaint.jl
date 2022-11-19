using HDF5
using Healpix
using Random
using Random: MersenneTwister


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
   Threads.@threads for (i1, i2) in chunk(num, chunksize)
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
Fill in basic halo properties.
"""
function get_basic_halo_properties(halo_pos::Array{T,2}, model::AbstractForegroundModel,
                                   cosmo::Cosmology.FlatLCDM{T}, res::Resolution) where T
    N_halos = size(halo_pos, 2)
    hp_ind = Array{Int64}(undef, N_halos)  # healpix index of halo
    redshift = Array{T}(undef, N_halos)
    dist = Array{T}(undef, N_halos)

    r2z = build_r2z_interpolator(
        model.min_redshift, model.max_redshift, cosmo)
    Threads.@threads for i in 1:N_halos
        dist[i] = sqrt(halo_pos[1,i]^2 + halo_pos[2,i]^2 + halo_pos[3,i]^2)
        redshift[i] = r2z(dist[i])
        hp_ind[i] = Healpix.vec2pixRing(res, halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
    end

    return dist, redshift, hp_ind
end

"""
Compute angles of halos
"""
function get_angles(halo_pos::Array{T,2}) where T
    N_halos = size(halo_pos, 2)
    θ = Array{T}(undef, N_halos)
    ϕ = Array{T}(undef, N_halos)

    Threads.@threads for i in 1:N_halos
        θ[i], ϕ[i] = Healpix.vec2ang(halo_pos[1,i], halo_pos[2,i], halo_pos[3,i])
    end

    return θ, ϕ
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


function catalog2map!(m::HealpixMap{T,RingOrder}, flux, theta, phi) where T
    res = m.resolution
    pixel_array = m.pixels
    N_halo = length(flux)

    # try to prevent thread issues by sorting by theta
    perm = sortperm(theta, rev=true, alg=ThreadsX.MergeSort)
    Threads.@threads for i_perm in 1:N_halo
        i_halo = perm[i_perm]
        hp_ind = Healpix.ang2pixRing(res, theta[i_halo], phi[i_halo])
        pixel_array[hp_ind] += flux[i_halo]
    end

    # divide by healpix pixel size
    per_pixel_steradian = 1 / nside2pixarea(res.nside)
    pixel_array .*= per_pixel_steradian
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
    Threads.@threads for i in 1:npix
        posmap.pixels[i] = pix2vecRing(res, i)
    end
    return posmap
end
vectorhealpixmap(nside::Int) = vectorhealpixmap(Float64, nside)


"""
Computes a real-space beam interpolator and a maximum
"""
function realspacegaussbeam(::Type{T}, θ_FWHM::Ti; rtol=1e-24, N_θ::Int=2000) where {T,Ti}
    Nlmax = ceil(Int, log2(8π / θ_FWHM))
    lmax = 2^Nlmax

    b_l = gaussbeam(θ_FWHM, lmax)
    θs = LinRange(zero(θ_FWHM), 5θ_FWHM, N_θ)
    b_θ = XGPaint.bl2beam(b_l, θs)
    atol = b_θ[begin] * rtol
    i_max = findfirst(<(atol), b_θ)

    θs = convert(LinRange{T, Int}, θs[begin:i_max])
    beam_real_interp = cubic_spline_interpolation(
        θs, T.(b_θ[begin:i_max]))
    return beam_real_interp, θs
end


struct BeamPaintingWorkspace{T, BI}
    ringinfo::RingInfo
    disc_buffer::Vector{Int}
    θmax::T
    posmap::HealpixMap{Tuple{T,T,T}, RingOrder}
    beam_real_interp::BI
end

function BeamPaintingWorkspace(nside::Int, θmax::T, beam_interp::BI) where {T, BI}
    ringinfo = RingInfo(0, 0, 0, 0.0, true)
    disc_buffer = Int[]
    approx_size = ceil(Int, 1.1 * π * θmax^2 / (nside2pixarea(nside)))  # 1.1 is fudge factor
    sizehint!(disc_buffer, approx_size)
    posmap = vectorhealpixmap(T, nside)
    return BeamPaintingWorkspace{T, BI}(ringinfo, disc_buffer, θmax, posmap, beam_interp)
end

function realspacebeampaint!(hp_map, w::BeamPaintingWorkspace, flux, θ₀, ϕ₀)
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    XGPaint.queryDiscRing!(w.disc_buffer, w.ringinfo, hp_map.resolution, θ₀, ϕ₀, w.θmax)

    for ir in w.disc_buffer
        x₁, y₁, z₁ = w.posmap[ir]
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(1 - d² / 2)

        hp_map.pixels[ir] += flux * w.beam_real_interp(θ)
    end
end


export read_halo_catalog_hdf5, chunk, generate_subhalo_offsets

