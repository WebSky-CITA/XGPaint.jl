"""
Map utilities for rectangular pixels.
"""

using WCS  # wraps the WCSLIB C library
using FFTW

"""
Map type, contains an AbstractArray and a WCS object, but behaves like the
AbstractArray it contains for array operations.

It only implements the subset of Base.Array operations which are common on maps.
You should work with the data directly using `enmap_instance.data` if you need
additional Array functions.
"""
struct Enmap{T,N,AA<:AbstractArray} <: AbstractArray{T,N}
    data::AA  # some kind of abstract array
    wcs::WCSTransform  # WCS object from WCS.jl
end

# getters
Base.parent(A::Enmap) = A.data
get_wcs(A::Enmap) = A.wcs

# constructor
Enmap(A::AbstractArray{T,N}, wcs::WCSTransform) where {T,N} =
    Enmap{T,N,typeof(A)}(A, wcs)

# array details
parenttype(::Type{Enmap{T,N,AA}}) where {T,N,AA} = AA
parenttype(A::Enmap) = parenttype(typeof(A))
Base.eachindex(::IndexCartesian, A::Enmap) = CartesianIndices(axes(A))
Base.IndexStyle(::Type{EM}) where {EM<:Enmap} = IndexStyle(parenttype(EM))
Base.size(A::Enmap) = size(parent(A))
Base.size(A::Enmap, d) = size(parent(A), d)

# create a similar map by making a deep copy of the WCS
function Base.similar(A::Enmap, ::Type{T}, dims::Dims) where T
    B = Enmap(similar(parent(A), T, dims), deepcopy(get_wcs(A)))
end

# send all the get and set to the data inside the enmap
@inline function Base.getindex(A::Enmap{T,N,AA}, I::Vararg{Int,N}) where {T,N,AA}
    parent(A)[I...]
end
@inline function Base.getindex(A::Enmap{T,N,AA}, i::Int) where {T,N,AA}
    parent(A)[i]
end
@inline function Base.setindex!(A::Enmap{T,N,AA}, val, I::Vararg{Int,N}) where {T,N,AA}
    parent(A)[I...] = val
end
@inline function Base.setindex!(A::Enmap{T,N,AA}, val, i::Int) where {T,N,AA}
    parent(A)[i] = val
end


# area functions
function extent_intermediate(shape::Tuple, wcs::WCSTransform)
    return abs.(deg2rad.(shape .* wcs.cdelt))
end

function area(shape::Tuple, wcs::WCSTransform)
    return prod(extent_intermediate(shape, wcs))
end

function pixsize(shape::Tuple, wcs::WCSTransform)
    return area(shape, wcs) / prod(shape)
end

function pixsize(m::Enmap)
    return pixsize(size(m), m.wcs)
end

function get_fft_norm(x::Enmap)
    return sqrt(pixsize(x) / prod(size(x.data)))
end

# normalizations made to agree with enmap
function enfft(x::Enmap{T,N,AA}) where {T,N,AA}
    return Enmap(FFTW.fft(x.data) .* get_fft_norm(x),
        deepcopy(x.wcs))
end

function enifft(x::Enmap{T,N,AA}) where {T,N,AA}
    return Enmap(FFTW.ifft(x.data) .* (1.0 / get_fft_norm(x)), deepcopy(x.wcs))
end

"""
Physically normalized enmap FFT
"""
function enfft!(x::Enmap{T,N,AA}) where {T,N,AA}
    FFTW.fft!(x.data)
    x .*= get_fft_norm(x)
end

"""
Physically normalized enmap inverse FFT
"""
function enifft!(x::Enmap{T,N,AA}) where {T,N,AA}
    FFTW.ifft!(x.data)
    x .*= (1.0/get_fft_norm(x))
end
