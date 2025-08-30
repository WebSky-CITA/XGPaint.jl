using XGPaint
using XGPaint: RingWorkspace
using Healpix
using Printf
using Test

"""
Brute-force reference implementation that loops over all pixels
"""
function profile_paint_bruteforce!(m::HealpixMap{T, RingOrder}, model, Mh, z, α₀, δ₀, θmax, normalization=1) where T
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

"""
Test the RingWorkspace implementation against brute force
"""
function test_ringworkspace_vs_bruteforce()
    @printf "Testing RingWorkspace vs brute force implementation...\n"
    
    # Setup
    nside = 64
    res = Healpix.Resolution(nside)
    workspace = RingWorkspace(res)
    
    # Create test maps
    map_ringworkspace = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
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
    
    all_passed = true
    
    for (i, test_case) in enumerate(test_cases)
        @printf "Test case %d: α₀=%.3f, δ₀=%.3f, θmax=%.3f\n" i test_case.α₀ test_case.δ₀ test_case.θmax
        
        # Reset maps
        fill!(map_ringworkspace.pixels, 0.0)
        fill!(map_bruteforce.pixels, 0.0)
        
        # Paint with both methods
        XGPaint.profile_paint_generic!(map_ringworkspace, workspace, model, 
                                     test_case.Mh, test_case.z, test_case.α₀, test_case.δ₀, test_case.θmax)
        
        profile_paint_bruteforce!(map_bruteforce, model, 
                                test_case.Mh, test_case.z, test_case.α₀, test_case.δ₀, test_case.θmax)
        
        # Compare results
        diff = map_ringworkspace.pixels .- map_bruteforce.pixels
        max_diff = maximum(abs.(diff))
        nonzero_ringworkspace = count(x -> x != 0, map_ringworkspace.pixels)
        nonzero_bruteforce = count(x -> x != 0, map_bruteforce.pixels)
        
        @printf "  Nonzero pixels: RingWorkspace=%d, BruteForce=%d\n" nonzero_ringworkspace nonzero_bruteforce
        @printf "  Max difference: %.2e\n" max_diff
        
        # Check if results are identical (should be exactly equal for constant model)
        if max_diff > 1e-14 || nonzero_ringworkspace != nonzero_bruteforce
            @printf "  ❌ FAILED: Results differ!\n"
            
            # Print some details about the differences
            nonzero_diff_count = count(x -> abs(x) > 1e-14, diff)
            @printf "  Pixels with significant differences: %d\n" nonzero_diff_count
            
            if nonzero_diff_count > 0 && nonzero_diff_count < 20
                for (idx, d) in enumerate(diff)
                    if abs(d) > 1e-14
                        @printf "    Pixel %d: RingWorkspace=%.6f, BruteForce=%.6f, diff=%.2e\n" idx map_ringworkspace.pixels[idx] map_bruteforce.pixels[idx] d
                    end
                end
            end
            
            all_passed = false
        else
            @printf "  ✅ PASSED\n"
        end
        
        @printf "\n"
    end
    
    return all_passed
end

"""
Test edge cases and boundary conditions
"""
function test_edge_cases()
    @printf "Testing edge cases...\n"
    
    nside = 64
    res = Healpix.Resolution(nside)
    workspace = RingWorkspace(res)
    
    map_test = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
    model = TestModel(1.0)
    
    edge_cases = [
        (α₀=0.0, δ₀=1.5708-0.01, θmax=0.1, name="Close to north pole"),
        (α₀=0.0, δ₀=-1.5708+0.01, θmax=0.1, name="Close to south pole"),
        (α₀=0.0, δ₀=0.0, θmax=1e-6, name="Tiny disc"),
        (α₀=0.0, δ₀=0.0, θmax=3.0, name="Large disc"),
        (α₀=6.28318-1e-6, δ₀=0.0, θmax=0.1, name="Right at RA wrap"),
    ]
    
    all_passed = true
    
    for (i, test_case) in enumerate(edge_cases)
        @printf "Edge case %d: %s\n" i test_case.name
        
        fill!(map_test.pixels, 0.0)
        
        try
            XGPaint.profile_paint_generic!(map_test, workspace, model, 
                                         1e14, 0.5, test_case.α₀, test_case.δ₀, test_case.θmax)
            
            nonzero_count = count(x -> x != 0, map_test.pixels)
            @printf "  Nonzero pixels: %d\n" nonzero_count
            @printf "  ✅ PASSED (no errors)\n"
        catch e
            @printf "  ❌ FAILED with error: %s\n" e
            all_passed = false
        end
        
        @printf "\n"
    end
    
    return all_passed
end

"""
Performance comparison test
"""
function test_performance()
    @printf "Performance comparison...\n"
    
    nside = 256  # Larger map for meaningful timing
    res = Healpix.Resolution(nside)
    workspace = RingWorkspace(res)
    
    map_ringworkspace = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
    map_bruteforce = HealpixMap{Float64, RingOrder}(zeros(Float64, nside2npix(nside)))
    
    model = TestModel(1.0)
    
    # Test case
    α₀, δ₀, θmax = 1.5708, 0.7854, 0.1
    Mh, z = 1e14, 0.5
    
    # Warm up
    XGPaint.profile_paint_generic!(map_ringworkspace, workspace, model, Mh, z, α₀, δ₀, θmax)
    profile_paint_bruteforce!(map_bruteforce, model, Mh, z, α₀, δ₀, θmax)
    
    # Time RingWorkspace
    fill!(map_ringworkspace.pixels, 0.0)
    time_ringworkspace = @elapsed XGPaint.profile_paint_generic!(map_ringworkspace, workspace, model, Mh, z, α₀, δ₀, θmax)
    
    # Time brute force
    fill!(map_bruteforce.pixels, 0.0)
    time_bruteforce = @elapsed profile_paint_bruteforce!(map_bruteforce, model, Mh, z, α₀, δ₀, θmax)
    
    speedup = time_bruteforce / time_ringworkspace
    
    @printf "RingWorkspace time: %.6f seconds\n" time_ringworkspace
    @printf "Brute force time:   %.6f seconds\n" time_bruteforce
    @printf "Speedup:            %.1fx\n" speedup
    
    return speedup > 1.0  # Should be faster
end

"""
Run all tests
"""
function run_all_tests()
    @printf "%s\n" ("="^60)
    @printf "Testing XGPaint RingWorkspace Implementation\n"
    @printf "%s\n\n" ("="^60)
    
    test1_passed = test_ringworkspace_vs_bruteforce()
    @printf "\n%s\n\n" ("="^60)
    
    test2_passed = test_edge_cases()
    @printf "\n%s\n\n" ("="^60)
    
    test3_passed = test_performance()
    @printf "\n%s\n\n" ("="^60)
    
    if test1_passed && test2_passed && test3_passed
        @printf "🎉 ALL TESTS PASSED! 🎉\n"
        return true
    else
        @printf "❌ SOME TESTS FAILED\n"
        @printf "Correctness: %s\n" (test1_passed ? "✅" : "❌")
        @printf "Edge cases:  %s\n" (test2_passed ? "✅" : "❌") 
        @printf "Performance: %s\n" (test3_passed ? "✅" : "❌")
        return false
    end
end

# Run the tests
if abspath(PROGRAM_FILE) == @__FILE__
    run_all_tests()
end
