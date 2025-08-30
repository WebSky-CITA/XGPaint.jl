
# Ring workspace for fast disc queries; this is designed for CPU for now
struct HealpixRingProfileWorkspace{T<:AbstractFloat}
    res::Healpix.Resolution
    ring_thetas::Vector{T}     # Colatitude of each ring
    ring_shifts::Vector{T}     # Shift offset for each ring (0.0 or 0.5)
    ring_num_pixels::Vector{Int} # Number of pixels in each ring
    ring_first_pixels::Vector{Int} # First pixel index for each ring
end

function HealpixRingProfileWorkspace{T}(res::Healpix.Resolution) where T<:AbstractFloat
    num_rings = res.nsideTimesFour - 1
    
    ring_thetas = Vector{T}(undef, num_rings)
    ring_shifts = Vector{T}(undef, num_rings)
    ring_num_pixels = Vector{Int}(undef, num_rings)
    ring_first_pixels = Vector{Int}(undef, num_rings)
    
    for ring_idx in 1:num_rings
        ringinfo = Healpix.getringinfo(res, ring_idx)
        ring_thetas[ring_idx] = T(ringinfo.colatitude_rad)
        ring_shifts[ring_idx] = ringinfo.shifted ? T(0.5) : T(0.0)
        ring_num_pixels[ring_idx] = ringinfo.numOfPixels
        ring_first_pixels[ring_idx] = ringinfo.firstPixIdx
    end
    
    HealpixRingProfileWorkspace{T}(res, ring_thetas, ring_shifts, ring_num_pixels, ring_first_pixels)
end

HealpixRingProfileWorkspace(res::Healpix.Resolution) = HealpixRingProfileWorkspace{Float64}(res)

function get_relevant_rings(res::Healpix.Resolution, θ_center, radius)
    z_top = cos(max(0, θ_center - radius))
    z_bottom = cos(min(π, θ_center + radius))
    
    ring_top = Healpix.ringAbove(res, z_top)
    ring_bottom = Healpix.ringAbove(res, z_bottom) + 1
    
    return max(1, ring_top), min(res.nsideTimesFour - 1, ring_bottom)
end

"""
    get_ring_disc_ranges(workspace, ring_idx, θ_center, ϕ_center, radius)

Returns (range1, range2) of pixel indices on a ring that lie within the disc.

The spherical distance between points (θ₁,φ₁) and (θ₂,φ₂) is given by:
    cos(d) = cos(θ₁)cos(θ₂) + sin(θ₁)sin(θ₂)cos(φ₁ - φ₂)

For a point to be included in the disc: d ≤ radius.

Returns two ranges to handle φ wraparound at the 0/2π boundary:
- range1: Main contiguous range of pixels  
- range2: Additional range when disc crosses φ=0/2π seam (empty if no crossing)
"""
function get_ring_disc_ranges(workspace::HealpixRingProfileWorkspace{T}, ring_idx::Int, 
                              center_theta::T, center_phi::T, radius::T) where T
    # Validate that phi is in the valid range [0, 2π)
    @assert 0 <= center_phi < T(2π) "center_phi must be in the range [0, 2π), got $(center_phi)"
    
    nr = workspace.ring_num_pixels[ring_idx]
    ring_theta = workspace.ring_thetas[ring_idx]
    shift = workspace.ring_shifts[ring_idx]
    
    # Polar cap regions require special handling
    is_north_cap = ring_idx ≤ workspace.res.nside
    is_south_cap = ring_idx ≥ 3 * workspace.res.nside
    
    if is_north_cap || is_south_cap
        # Exact pole centers: simple distance check
        if abs(center_theta) < T(1e-10)
            return ring_theta ≤ radius ? (1:nr, 1:0) : (1:0, 1:0)
        end
        if abs(center_theta - T(π)) < T(1e-10)
            return (T(π) - ring_theta) ≤ radius ? (1:nr, 1:0) : (1:0, 1:0)
        end
        
        # Solve spherical triangle for phi half-width
        cos_center, sin_center = cos(center_theta), sin(center_theta)
        cos_ring, sin_ring = cos(ring_theta), sin(ring_theta)
        cos_radius = cos(radius)
        
        denominator = sin_ring * sin_center
        if abs(denominator) < T(1e-12)
            condition = abs(cos_radius - cos_ring * cos_center) < T(1e-12)
            return condition ? (1:nr, 1:0) : (1:0, 1:0)
        end
        
        cos_delta_phi = (cos_radius - cos_ring * cos_center) / denominator
        if cos_delta_phi > T(1); return (1:0, 1:0); end
        if cos_delta_phi < T(-1); return (1:nr, 1:0); end
        
        delta_phi = acos(cos_delta_phi)
        phi_min, phi_max = center_phi - delta_phi, center_phi + delta_phi
        
        # Convert phi to pixel indices
        ip_lo = floor(Int, nr / T(2π) * phi_min - shift) + 1
        ip_hi = floor(Int, nr / T(2π) * phi_max - shift)
        if ip_lo > ip_hi; return (1:0, 1:0); end
        
        # Handle phi wraparound at 0/2π boundary
        if ip_hi >= nr; ip_lo -= nr; ip_hi -= nr; end
        return ip_lo < 0 ? (1:(ip_hi + 1), (ip_lo + nr + 1):nr) : ((ip_lo + 1):(ip_hi + 1), 1:0)
    else
        # Equatorial region: standard HEALPix algorithm
        z, z0 = cos(ring_theta), cos(center_theta)
        xa = T(1) / sqrt((T(1) - z0) * (T(1) + z0))
        x = (cos(radius) - z * z0) * xa
        ysq = T(1) - z^2 - x^2
        
        dphi = ysq < T(0) ? T(0) : atan(sqrt(ysq), x)
        if dphi ≤ T(0); return (1:0, 1:0); end
        
        ip_lo = floor(Int, nr / T(2π) * (center_phi - dphi) - shift) + 1
        ip_hi = floor(Int, nr / T(2π) * (center_phi + dphi) - shift)
        
        if ip_lo > ip_hi; return (1:0, 1:0); end
        if ip_hi >= nr; ip_lo -= nr; ip_hi -= nr; end
        
        return ip_lo < 0 ? (1:(ip_hi + 1), (ip_lo + nr + 1):nr) : ((ip_lo + 1):(ip_hi + 1), 1:0)
    end
end
