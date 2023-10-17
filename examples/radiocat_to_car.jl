using HDF5
using Pixell
using Healpix
using XGPaint

cat_dir = "/global/project/projectdirs/sobs/www/users/Radio_WebSky/matched_catalogs_2"
# map_dir = "/global/project/projectdirs/sobs/www/users/Radio_WebSky/matched_maps/"
map_dir = "/global/cscratch1/sd/xzackli/radio_maps_will/"


shape, wcs = fullsky_geometry(deg2rad(0.5 / 60))
m = Enmap(zeros(shape), wcs)


## 
function enmap_py2jl(m_py)
    m0 = pyconvert(Array{Float64, 2}, m_py)
    permutedims(m0, (2,1))
end

using ThreadsX
using PythonCall
enmap = pyimport("pixell.enmap")

shape_py, wcs_py = enmap.fullsky_geometry(res=deg2rad(0.5 / 60))
pixsizes_py = enmap.pixsizemap(shape_py, wcs_py)
pixsizes = enmap_py2jl(pixsizes_py)

##
function catalog2map!(m::Enmap{T}, flux, theta, phi, pixsizes) where T
    N_halo = length(flux)
    pixel_array = parent(m)
    fill!(pixel_array, zero(T))
    α = phi
    δ = π/2 .- theta

    ipix, jpix = sky2pix(m, α, δ)

    # try to prevent thread issues by sorting by theta
    perm = sortperm(theta, rev=true, alg=ThreadsX.MergeSort)
    Threads.@threads :static for i_perm in 1:N_halo
        i_halo = perm[i_perm]
        i = round(Int, ipix[i_halo])
        i = mod(i - 1, size(m, 1)) + 1
        j = round(Int, jpix[i_halo])
        pixel_array[i, j] += flux[i_halo]
    end
    
    pixel_array ./= pixsizes
    
    top_pole = sum(pixel_array[:,end])
    bot_pole = sum(pixel_array[:,1])

    pixel_array[:,end] .= top_pole
    pixel_array[:,1] .= bot_pole
    return m
end

##

# filename = "/global/project/projectdirs/sobs/www/users/Radio_WebSky/matched_catalogs_2/catalog_90.2.h5"
# sources = convert(Dict{String, Vector{Float64}}, read(h5open(filename, "r")))
# catalog2map!(m, sources["flux"], sources["theta"], sources["phi"], pixsizes)
# write_map(map_path, m)

##
write_map("test.fits", m)

##
for (root, dirs, files) in walkdir(cat_dir)
    # println("Directories in $root")
    # for dir in dirs
    #     println(joinpath(root, dir)) # path to directories
    # end
    println("Files in $root")
    for file in files
        freq = replace(file, ".h5" => "", "catalog_" => "")
        filename = joinpath(root, file)
        println(filename)
        sources = convert(Dict{String, Vector{Float64}}, read(h5open(filename, "r")))
        catalog2map!(m, sources["flux"], sources["theta"], sources["phi"], pixsizes)
        map_path = joinpath(map_dir, "map_car_0.5arcmin_f$(freq).fits")
        write_map(map_path, m)
    end
end

##

wy, wx = enmap.calc_window(shape_py)




# function catalog2map!(m::HealpixMap{T,RingOrder}, flux, theta, phi) where T
#     res = m.resolution
#     pixel_array = m.pixels
#     N_halo = length(flux)

#     # try to prevent thread issues by sorting by theta
#     perm = sortperm(theta, rev=true, alg=ThreadsX.MergeSort)
#     Threads.@threads :static for i_perm in 1:N_halo
#         i_halo = perm[i_perm]
#         hp_ind = Healpix.ang2pixRing(res, theta[i_halo], phi[i_halo])
#         pixel_array[hp_ind] += flux[i_halo]
#     end

#     # divide by healpix pixel size
#     per_pixel_steradian = 1 / nside2pixarea(res.nside)
#     pixel_array .*= per_pixel_steradian
# end
