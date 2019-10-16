using XGPaint
using Cosmology
using Test

# relative background evolutions differ by 1e-3 between Julia and Python 2
rtol = 1e-3
cosmo = XGPaint.get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = XGPaint.CIBModel{Float32}()
r2z = XGPaint.build_r2z_interpolator(0.0f0, 4.5f0, cosmo)
hod_shang = XGPaint.build_shang_interpolator(log(1.0f13), log(1.0f15), model)
clnm2r = XGPaint.build_c_lnm2r_interp(nbin=30)

@testset "halo2fluxmap" begin
    # These tests come from computing examples in the Python xgpaint

    @test model.shang_Mpeak ≈ 10^12.3  # check default is viero as docs say

    # python: h2fm.utils.r2z(3407.0f0)
    @test isapprox(
        r2z(3407.0f0), 0.9996, rtol=rtol)

    # python: h2fm.utils.jiang_shmf(0.5e12, 1e13)
    @test XGPaint.jiang_shmf(0.5f12, 1.0f13, model) ≈ 2.5335820678141165

    # python: h2fm.hod.hod_shang(np.array([1e13, 1e14, 1e15]))
    @test hod_shang(log(1.0f13)) ≈ 9.78325852
    @test hod_shang(log(1.0f14)) ≈ 54.21666623
    @test hod_shang(log(1.0f15)) ≈ 336.55409152

    # python: h2fm.utils.mz2c(1e13, 4.0)
    @test isapprox(
        XGPaint.mz2c(1.0f13, 4.0f0, cosmo),
        2.262202312880881, rtol=rtol)

    # h2fm.utils.m2r(1e13)
    @test isapprox(
        XGPaint.m2r(1.0f13, cosmo),
        4.12328640425709, rtol=rtol)

    # h2fm.hod.clnm2r(np.array([1e-3, 25.0]), np.log(0.2))
    @test isapprox(
        clnm2r(1f-3 / 0.999f0, log(0.2f0)),
        0.44787404, rtol=rtol)
    @test isapprox(
        clnm2r(25.0f0 / 1.001f0, log(0.2f0)),
        0.085837, rtol=rtol)

    @test XGPaint.random_phi(Float32) <= 2π
    @test XGPaint.random_theta(Float32) <= π

    # h2fm.fluxmodel.sigma_cen(1e13)
    @test isapprox(XGPaint.sigma_cen(1.0f13, model),
        3218663770378.3066, rtol=rtol)

    # h2fm.fluxmodel.nu2theta(150e9, 0.5), h2fm.fluxmodel.nu2theta(1.5e6, 1.0)
    @test isapprox(XGPaint.nu2theta(150e9, 0.5, model),
        6.5705729824589e-16, rtol=rtol)
    @test isapprox(XGPaint.nu2theta(1.5e6, 1.0, model),
        5.945023942373837e-34, rtol=rtol)

    # h2fm.fluxmodel.integrand_L(30, 32)
    @test isapprox(XGPaint.integrand_L(30.0f0, 32.0f0, model),
        3929486334377.691, rtol=rtol)
end
