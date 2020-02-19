m = Map{Float64,RingOrder}(model.nside)

using BenchmarkTools
using HDF5
using JLD2

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
function testJLD2()
    A = m.pixels
    @save "test.jld2" A
end


##
@time testJLD2()
##
@time testHDF5()
##
@time Healpix.saveToFITS(m, "!test.fits", typechar="E")
##
@time Healpix.saveToFITS(m, "!test.fits", typechar="D")
##
