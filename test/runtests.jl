using XGPaint
using Cosmology
using Test

# relative background evolutions differ by 1e-3 between Julia and Python
rtol = 1e-3
cosmo = cosmology(h=0.7, OmegaM=0.25)
r2z = XGPaint.build_r2z_interpolator(0.0f0, 4.5f0, cosmo)
model = CIBModel{Float32}()
hod_shang = XGPaint.build_shang_interpolator(log(1.0f13), log(1.0f15), model)

@testset "halo2fluxmap" begin
    # These tests come from computing specific numbers from the Python xgpaint

    @test model.shang_Mpeak ≈ 10^12.3  # check default is viero as docs say

    # python: h2fm.utils.r2z(3407.0f0)
    @test isapprox(r2z(3407.0f0), 0.9996, rtol=rtol)

    # python: h2fm.utils.jiang_shmf(0.5e12, 1e13)
    @test XGPaint.jiang_shmf(0.5f12, 1.0f13, model) ≈ 2.5335820678141165

    # python: h2fm.hod.hod_shang(np.array([1e13, 1e14, 1e15]))
    @test hod_shang(log(1.0f13)) ≈ 9.78325852
    @test hod_shang(log(1.0f14)) ≈ 54.21666623
    @test hod_shang(log(1.0f15)) ≈ 336.55409152

end
