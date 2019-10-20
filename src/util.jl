using HDF5
using Random

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
Generate an array of random number generators, for each thread.
"""
function get_thread_RNG()
    r = let m = MersenneTwister(1)
            [m; accumulate(Future.randjump, fill(big(10)^20,
            Threads.nthreads()-1), init=m)]
        end;
    return r
end

"""
Generates a list of tuples which describe starting and ending chunk indices.
Useful for parallelizing an array operation.
"""
function chunk_list(arr_len, chunksize)
    return [(i,min(i + chunksize-1, arr_len))
        for i in range(1, arr_len, step=chunksize)]
end

export read_halo_catalog_hdf5, chunk_list
