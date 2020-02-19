m = Map{Float64,RingOrder}(model.nside)

using BenchmarkTools
using HDF5

@btime Healpix.saveToFITS(m, "!test.fits", typechar="E")

##

function testHDF5()
    h5open("test.h5", "w") do file
        write(file, "A", m.pixels)  # alternatively, say "@write file A"
    end
end

##
@btime testHDF5()


##
