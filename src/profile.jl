
"""
Perform an integral of the function `f(x,y)` over a rectangle defined by x_range and y_range,
using nested Gauss-Legendre quadrature. 
"""
function integrate_over_rect(::Type{T}, f, x_range, y_range, gl_grid, gl_weights) where T
    @assert length(gl_grid) == length(gl_weights)  # check we have sensible Gauss-Legendre xᵢ and wᵢ
    N = length(gl_grid)

    # integral is over (x₁, x₂), (y₁, y₂) instead of (-1,1), (-1,1), but grid and weights are in the latter
    x₁, x₂ = x_range
    y₁, y₂ = y_range
    m_x = (x₂ - x₁) / 2
    b_x = (x₁ + x₂) / 2
    m_y = (y₂ - y₁) / 2
    b_y = (y₁ + y₂) / 2
    change_of_variables = m_x * m_y
    
    result = zero(T)  # accumulate here
    @turbo for i in 1:N
        y_result = zero(T)
        for j in 1:N
            x = m_x * gl_grid[i] + b_x
            y = m_y * gl_grid[j] + b_y
            y_result += f(x,y) * gl_weights[j]
        end
        result += y_result * gl_weights[i]
    end

    return result * change_of_variables
end


abstract type ProfilePaintingMethod end
struct SamplePixelCenters <: ProfilePaintingMethod end


"""

"""
function paint_profile!(::SamplePixelCenters, m::Enmap, α, δ,
                        radial_profile, total_flux, buffer; safe=true)
    center = (α, δ)
    center_i, center_j = round.(Int64, sky2pix(m, α, δ; safe=safe))
    sky2pix(m, α, δ; safe=false)
    buf_center_i = size(buffer, 1) ÷ 2
    buf_center_j = size(buffer, 2) ÷ 2

    metric = SphericalAngle()  # using Haversine distance from Distances.jl, in radians

    # compute distance 
    for j in axes(buffer, 2)
        for i in axes(buffer, 1)
            # get pixel in sky coordinates, and then compute the arcelength on the sky
            map_i = center_i + i - buf_center_i
            map_j = center_j + j - buf_center_j
            pixloc = pix2sky(m, map_i, map_j)
            # Ω_p = pixsize(m, map_i, map_j)  # NEED TO IMPLEMENT pixsize
            arclength = metric(center, pixloc)
            buffer[i,j] = radial_profile(arclength) * Ω_p
        end
    end

    norm = total_flux / sum(buffer)

    # now add buffer to Enmap
    map_arr = parent(m)
    for j in axes(buffer, 2)
        for i in axes(buffer, 1)
            map_i, map_j = center_i + i - buf_center_i, center_j + j - buf_center_j
            if (1 ≤ map_i ≤ size(map_arr,1)) && (1 ≤ map_j ≤ size(map_arr,2))
                @inbounds map_arr[map_i, map_j] = buffer[i,j] * norm
            end
        end
    end

end
