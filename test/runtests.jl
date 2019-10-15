using XGPaint
using Cosmology
using Test

# background evolution is within 1e-7 between Julia and Python implementations
rtol = 1e-6
cosmo = cosmology(h=0.7, OmegaM=0.25)
r2z = XGPaint.build_r2z_interpolator(0.0f0, 4.5f0, cosmo)

@testset "halo2fluxmap" begin
    # These tests come from computing specific numbers from the Python xgpaint
    model = CIBModel{Float32}()
    @test model.shang_Mpeak ≈ 10^12.3  # check default is viero as docs say
    @test isapprox(r2z(3407.6934323399264f0), 1.0, rtol=rtol)

end
