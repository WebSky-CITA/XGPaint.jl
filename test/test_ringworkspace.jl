using XGPaint
using XGPaint: HealpixRingProfileWorkspace
using Healpix
using Test

"""
Brute-force reference implementation that loops over all pixels
"""
function profile_paint_bruteforce!(
    m::HealpixMap{T, RingOrder}, model, Mh, z, α₀, δ₀, θmax, normalization=1
) where T

    ϕ₀ = α₀  
    θ₀ = T(π)/2 - δ₀
    x₀, y₀, z₀ = ang2vec(θ₀, ϕ₀)
    θmin = XGPaint.compute_θmin(model)
    
    # Brute force: check every pixel in the map
    for global_pix in 1:length(m.pixels)
        # Get position of this pixel
        x₁, y₁, z₁ = pix2vecRing(m.resolution, global_pix)
        
        # Compute angular distance
        d² = (x₁ - x₀)^2 + (y₁ - y₀)^2 + (z₁ - z₀)^2
        θ = acos(clamp(1 - d² / 2, -one(T), one(T)))
        θ = max(θmin, θ)  # clamp to minimum θ
        
        # Add contribution to map
        m.pixels[global_pix] += ifelse(θ < θmax,
                                      normalization * model(θ, Mh, z),
                                      zero(T))
    end
end

"""
Simple test model that returns a constant value for testing
"""
struct TestModel{T} <: XGPaint.AbstractProfile{T}
    value::T
end

(m::TestModel)(θ, Mh, z) = m.value
XGPaint.compute_θmin(::TestModel{T}) where T = eps(T)

@testset "HealpixRingProfileWorkspace Profile Painting Tests" begin

    @testset "HealpixRingProfileWorkspace vs Brute Force" begin
        # Setup
        nside = 64
        res = Healpix.Resolution(nside)
        workspace = HealpixRingProfileWorkspace(res)
        
        # Create test maps
        map_HealpixRingProfileWorkspace = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
        map_bruteforce = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
        
        # Test model
        model = TestModel(1.0)
        
        # Test parameters
        test_cases = [
            (α₀=0.0, δ₀=0.0, θmax=0.1, Mh=1e14, z=0.5),           # North pole
            (α₀=3.14159, δ₀=0.0, θmax=0.1, Mh=1e14, z=0.5),       # South pole  
            (α₀=1.5708, δ₀=0.7854, θmax=0.2, Mh=1e14, z=0.3),     # Equatorial
            (α₀=4.7124, δ₀=-0.7854, θmax=0.15, Mh=1e13, z=0.7),   # Southern hemisphere
            (α₀=0.1, δ₀=0.05, θmax=0.05, Mh=5e13, z=1.0),         # Small disc
            (α₀=6.1832, δ₀=0.0, θmax=0.3, Mh=2e14, z=0.1),        # Near RA=0 wrap
        ]
        
        for (i, test_case) in enumerate(test_cases)
            # Reset maps
            fill!(map_HealpixRingProfileWorkspace.pixels, 0.0)
            fill!(map_bruteforce.pixels, 0.0)
            
            # Paint with both methods
            XGPaint.profile_paint_generic!(map_HealpixRingProfileWorkspace, workspace, model, 
                                         test_case.Mh, test_case.z, test_case.α₀, test_case.δ₀, test_case.θmax)
            
            profile_paint_bruteforce!(map_bruteforce, model, 
                                    test_case.Mh, test_case.z, test_case.α₀, test_case.δ₀, test_case.θmax)
            
            # Compare results - should be exactly equal for constant model
            @test map_HealpixRingProfileWorkspace.pixels ≈ map_bruteforce.pixels
            @test count(x -> x != 0, map_HealpixRingProfileWorkspace.pixels) == count(x -> x != 0, map_bruteforce.pixels)
        end
    end

    @testset "HealpixRingProfileWorkspace Edge Cases" begin
        nside = 64
        res = Healpix.Resolution(nside)
        workspace = HealpixRingProfileWorkspace(res)
        
        map_test = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
        model = TestModel(1.0)
        
        edge_cases = [
            (α₀=0.0, δ₀=1.5708-0.01, θmax=0.1, name="Close to north pole"),
            (α₀=0.0, δ₀=-1.5708+0.01, θmax=0.1, name="Close to south pole"),
            (α₀=0.0, δ₀=0.0, θmax=1e-6, name="Tiny disc"),
            (α₀=0.0, δ₀=0.0, θmax=3.0, name="Large disc"),
            (α₀=6.28318-1e-6, δ₀=0.0, θmax=0.1, name="Right at RA wrap"),
        ]
        
        for test_case in edge_cases
            fill!(map_test.pixels, 0.0)
            
            # Should not throw errors
            @test_nowarn XGPaint.profile_paint_generic!(map_test, workspace, model, 
                                                       1e14, 0.5, test_case.α₀, test_case.δ₀, test_case.θmax)
        end
    end

    @testset "HealpixRingProfileWorkspace vs Brute Force with Real Y Profile" begin
        # Setup with real y profile model
        nside = 32  # Smaller for faster testing with real profile
        res = Healpix.Resolution(nside)
        workspace = HealpixRingProfileWorkspace(res)
        
        # Create test maps
        map_HealpixRingProfileWorkspace = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
        map_bruteforce = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
        
        # Real y profile model
        y_model = XGPaint.Battaglia16ThermalSZProfile()
        
        # Test parameters - realistic values
        test_cases = [
            (α₀=0.0, δ₀=0.0, θmax=0.05, Mh=1e14, z=0.3),         # North pole, typical cluster
            (α₀=1.5708, δ₀=0.5236, θmax=0.08, Mh=5e13, z=0.7),   # Mid-sky, smaller cluster
            (α₀=3.14159, δ₀=-0.5236, θmax=0.03, Mh=2e14, z=0.1), # South, massive nearby cluster
        ]
        
        for (i, test_case) in enumerate(test_cases)
            # Reset maps
            fill!(map_HealpixRingProfileWorkspace.pixels, 0.0)
            fill!(map_bruteforce.pixels, 0.0)
            
            # Paint with both methods
            XGPaint.profile_paint_generic!(map_HealpixRingProfileWorkspace, workspace, y_model, 
                                         test_case.Mh, test_case.z, test_case.α₀, test_case.δ₀, test_case.θmax)
            
            profile_paint_bruteforce!(map_bruteforce, y_model, 
                                    test_case.Mh, test_case.z, test_case.α₀, test_case.δ₀, test_case.θmax)
            
            # Compare results - should be very close (allowing for floating point differences)
            @test map_HealpixRingProfileWorkspace.pixels ≈ map_bruteforce.pixels rtol=1e-12
            
            # Check that we actually painted something non-zero
            @test count(x -> x != 0, map_HealpixRingProfileWorkspace.pixels) > 0
            @test count(x -> x != 0, map_bruteforce.pixels) > 0
            
            # Check same number of non-zero pixels
            @test count(x -> x != 0, map_HealpixRingProfileWorkspace.pixels) == count(x -> x != 0, map_bruteforce.pixels)
        end
    end

end  # HealpixRingProfileWorkspace Profile Painting Tests
