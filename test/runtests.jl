using XGPaint
using Test

@test f(1) == 2

@testset "XGPaint.jl" begin
    # Write your own tests here.
    @test f(1) == 2
end
