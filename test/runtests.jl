using XGPaint
using Cosmology
using Test
using HDF5
using Pixell

using XGPaint: M_sun, get_params, generalized_nfw, rho_2d, ne2d
using Unitful, UnitfulAstro
using PhysicalConstants.CODATA2018

# relative background evolutions differ by 1e-3 between Julia and Python 2
rtol = 1e-3
cosmo = XGPaint.get_cosmology(h=0.7f0, OmegaM=0.25f0)

@testset "cib" begin
    # These tests come from computing examples in the Python xgpaint
    model = XGPaint.CIB_Planck2013{Float32}()
    r2z = XGPaint.build_r2z_interpolator(0.0f0, 4.5f0, cosmo)
    hod_shang = XGPaint.build_shang_interpolator(log(1.0f13), log(1.0f15), model)
    clnm2r = XGPaint.build_c_lnm2r_interpolator(nbin=30)
    #abandon sigma_sat for now
    #sigma_sat = XGPaint.build_sigma_sat_ln_interpolator(log(4f15), model)
    muofn = XGPaint.build_muofn_interpolator(model)

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

    # python: h2fm.utils.m2r(1e13)
    @test isapprox(
        XGPaint.m2r(1.0f13, cosmo),
        4.12328640425709, rtol=rtol)

    # python: h2fm.hod.clnm2r(np.array([1e-3, 25.0]), np.log(0.2))
    @test isapprox(
        clnm2r(1f-3 / 0.999f0, log(0.2f0)),
        0.44787404, rtol=rtol)
    @test isapprox(
        clnm2r(25.0f0 / 1.001f0, log(0.2f0)),
        0.085837, rtol=rtol)

    # python: h2fm.fluxmodel.sigma_cen(1e13)
    @test isapprox(XGPaint.sigma_cen(1.0f13, model),
        3218663770378.3066, rtol=rtol)

    # python: h2fm.fluxmodel.nu2theta(150e9, 0.5)
    @test isapprox(XGPaint.nu2theta(150f9, 0.5f0, model),
        6.5705729824589e-16/2, rtol=rtol)
    # python: h2fm.fluxmodel.nu2theta(1.5e6, 1.0)
    @test isapprox(XGPaint.nu2theta(1.5e6, 1.0, model),
        5.945023942373837e-34/2, rtol=rtol)

    #abandon sigma_sat for now
    # python: h2fm.fluxmodel.integrand_L(30, 32)
    #@test isapprox(
    #    XGPaint.integrand_L(30.0f0, 32.0f0, model),
    #    3929486334377.691, rtol=rtol)

    # python: h2fm.fluxmodel.l2f(1.0, 1.0e-3, 1.0e-3, 1.0e-3)
    @test isapprox(
        XGPaint.l2f( 1.0f0, sqrt(3.0f-6), r2z(sqrt(3.0f-6))) * 4π,
        333198.42738065525, rtol=rtol)
    #abandon sigma_sat for now
    # python: h2fm.fluxmodel.sigma_sat(np.array([1e13, 4e15]))
    #@test sigma_sat(log(4f15)) ≈ 3.50530158e+14
    #@test sigma_sat(log(1f13)) ≈ 2.48648819e+12

    @test muofn(500.0f0) ≈ 6.14975653e-05
    @test muofn(1.0f0) ≈ 0.11765558

    @test XGPaint.z_evo(0.0f0, model) ≈ 1.0f0
end

@testset "radio" begin
    radio_model = Radio_Sehgal2009{Float32}()
    @test XGPaint.FR_I_redshift_evolution(0.0f0, radio_model) ≈ 1.0f0
    @test XGPaint.FR_II_redshift_evolution(0.0f0, radio_model) ≈ 1.0f0
    @test (XGPaint.FR_I_redshift_evolution(radio_model.I_z_p, radio_model) ≈
        XGPaint.FR_I_redshift_evolution(radio_model.I_z_p+1.0f0, radio_model))
    @test XGPaint.FR_II_redshift_evolution(1.3f0, radio_model) > 190

end

@testset "sanity" begin
    @test XGPaint.l2f(1.0, 1.0, 0.0) ≈ 1/4π
    @test XGPaint.random_phi(Float32) <= 2π
    @test XGPaint.random_theta(Float32) <= π
    @test XGPaint.random_phi(Float64) <= 2π
    @test XGPaint.random_theta(Float64) <= π

    # file reading test
    testData = rand( Float32, 4, 100 )
    h5open("test_file_writing.h5", "w") do file
        write(file, "halos", testData)  # alternatively, say "@write file A"
    end
    read_pos, read_mass = XGPaint.read_halo_catalog_hdf5("test_file_writing.h5")
    @test all( read_pos .≈ testData[1:3,:] )
    @test all( read_mass .≈ testData[4,:] )
    rm("test_file_writing.h5")

    @test all(XGPaint.chunk(10, 3) == [(1,3), (4,6), (7,9), (10,10)])
    @test all(XGPaint.chunk(10, 4) == [(1,4), (5,8), (9,10)])
    @test all(XGPaint.chunk(10, 10) == [(1,10)])

end


include("test_query.jl")

@testset "tsz" begin
    ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
    ra, dec, redshift, halo_mass = sort_halo_catalog(ra, dec, redshift, halo_mass);
    y_model_interp = XGPaint.load_precomputed_battaglia()
    box = [4.5   -4.5;           # RA
       -3     3] * Pixell.degree  # DEC
    shape, wcs = geometry(Pixell.CarClenshawCurtis{Float64}, box, 0.5 * Pixell.arcminute)
    m = Enmap(zeros(shape), wcs)
    workspace = profileworkspace(shape, wcs)
    @time paint!(m, workspace, y_model_interp, halo_mass, redshift, ra, dec)

    ymap_ref = XGPaint.load_example_tsz_map()
    @test sum(abs, m.data .- ymap_ref) ≈ 0 atol=1e-7
end

@testset "tau_profile" begin
    # test the 3D tau profile
    p = XGPaint.BattagliaTauProfilePhysical(Omega_c=0.267, Omega_b=0.0493,  h=0.6712)

    # fits are in Msun/h, will change later
    par = get_params(p, (1e14M_sun / 0.6712), 0.5)
    @test par.P₀ * generalized_nfw(0.5, par.xc, par.α, par.β, par.γ) ≈ 231.94578059850758

    par = get_params(p, (1e14M_sun / 0.6712), 0.8)
    @test par.P₀ * generalized_nfw(0.5, par.xc, par.α, par.β, par.γ) ≈ 228.24947739841136

    par = get_params(p, (1e15M_sun / 0.6712), 0.8)
    @test par.P₀ * generalized_nfw(0.2, par.xc, par.α, par.β, par.γ) ≈ 1447.2222096644925

    par = get_params(p, (1e15M_sun / p.cosmo.h), 0.8)
    @test par.P₀ * XGPaint._nfw_profile_los_quadrature(0.5,  par.xc, par.α, par.β, par.γ) ≈ 320.36848661635133


    zz = 0.5

    Mnew = 86349927525539.45 * (M_sun / p.cosmo.h)
    ta = rho_2d(p, 1u"Mpc" / (1+zz) / p.cosmo.h, Mnew, zz) / (M_sun / 1u"Mpc^2") * p.cosmo.h + 0
    @test abs(1 - ta / 10.149657232478662e12) < 1e-4

    tb = rho_2d(p, 2u"Mpc" / (1+zz) / p.cosmo.h, Mnew, zz) / (M_sun / 1u"Mpc^2") * p.cosmo.h + 0
    @test abs(1 - tb / 2.446742995577804e12) < 1e-4

    zz = 2.5
    Mnew = 93218298413772.23 * (M_sun / p.cosmo.h)  # Mcrit
    tc = rho_2d(p, 3u"Mpc" / (1+zz) / p.cosmo.h, Mnew, zz) / (M_sun / 1u"Mpc^2") * p.cosmo.h + 0
    @test abs(1 - tc / 0.9123565449059163e12) < 1e-4

    zz = 2.5
    Mnew = 93218298413772.23 * (M_sun / p.cosmo.h)  # Mcrit
    tc = ne2d(p, 3u"Mpc" / (1+zz) / p.cosmo.h, Mnew, zz) / (p.cosmo.h / 1u"Mpc")^2 + 0
    @test abs(1 - tc / 1.4217533501414173e+69) < 1e-3

    zz = 2.5
    Mnew = 93218298413772.23 * (XGPaint.M_sun / p.cosmo.h)  # Mcrit
    tc = XGPaint.tau(p, 3u"Mpc" / (1+zz) / p.cosmo.h, Mnew, zz) + 0
    @test abs(1 - tc / 4.475127577749756e-05) < 1e-3

end
