using XGPaint
using Test

@testset "XGPaint.jl" begin
    # Write your own tests here.
    # @test f(1) == 2
    model = CIBModel{Float32}()
    @test model.min_mass ≈ 1.0f12  # check default is viero as docs say
end
