##
using XGPaint, Random, PyPlot, BenchmarkTools
#
# function threaded_rand!(random_number_generators, arr::Array{T,1};
#       chunksize=4096) where T
#
#    num = size(arr,1)
#    Threads.@threads for (i1, i2) in chunk(num, chunksize)
#       @views rand!(random_number_generators[Threads.threadid()], arr[i1:i2])
#    end
# end

rng_array = get_thread_RNG()

##

v = -1 .* ones(1000000000)
##
@btime XGPaint.threaded_rand!(rng_array, v; chunksize=4096)
##
@btime rand!(rng_array[1], v)
##
