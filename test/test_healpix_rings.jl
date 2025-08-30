using Test
using Healpix
using XGPaint

@testset "Ring-Slice HEALPix Disc Query Tests" begin
    
    @testset "HealpixRingProfileWorkspace Construction" begin
        nside = 32
        res = Healpix.Resolution(nside)
        workspace = XGPaint.HealpixRingProfileWorkspace(res)
        
        @test workspace.res === res
        @test length(workspace.ring_thetas) == res.nsideTimesFour - 1
        @test length(workspace.ring_shifts) == res.nsideTimesFour - 1
        @test length(workspace.ring_num_pixels) == res.nsideTimesFour - 1
        @test eltype(workspace.ring_thetas) == Float64
        @test eltype(workspace.ring_shifts) == Float64
        
        # Test type parameter
        workspace_f32 = XGPaint.HealpixRingProfileWorkspace{Float32}(res)
        @test eltype(workspace_f32.ring_thetas) == Float32
        @test eltype(workspace_f32.ring_shifts) == Float32
    end
    
    @testset "get_relevant_rings" begin
        nside = 32
        res = Healpix.Resolution(nside)
        
        # Test north pole
        ring_min, ring_max = XGPaint.get_relevant_rings(res, 0.0, deg2rad(10.0))
        @test ring_min == 1
        @test ring_max ≥ ring_min
        
        # Test equatorial region
        ring_min, ring_max = XGPaint.get_relevant_rings(res, π/2, deg2rad(5.0))
        @test 1 ≤ ring_min ≤ ring_max ≤ res.nsideTimesFour - 1
        
        # Test south pole
        ring_min, ring_max = XGPaint.get_relevant_rings(res, π, deg2rad(10.0))
        @test ring_min ≤ ring_max
        @test ring_max == res.nsideTimesFour - 1
    end
    
    @testset "get_ring_disc_ranges Basic Cases" begin
        nside = 32
        res = Healpix.Resolution(nside)
        workspace = XGPaint.HealpixRingProfileWorkspace(res)
        
        # Test empty disc (radius = 0)
        range1, range2 = XGPaint.get_ring_disc_ranges(workspace, 10, π/4, 0.0, 0.0)
        @test isempty(range1) || length(range1) ≤ 1
        @test isempty(range2)
        
        # Test large disc that includes significant portion of ring
        range1, range2 = XGPaint.get_ring_disc_ranges(workspace, 50, π/2, 0.0, deg2rad(45.0))
        expected_pixels = workspace.ring_num_pixels[50]
        # With a 45° radius, we should get a reasonable number of pixels
        @test length(range1) + length(range2) ≥ 1
    end
    
    @testset "Systematic Validation vs HEALPix Reference" begin
        nside = 32  # Smaller for faster testing
        res = Healpix.Resolution(nside)
        workspace = XGPaint.HealpixRingProfileWorkspace(res)
        
        function query_disc_rings!(pixel_list::Vector{Int}, workspace::XGPaint.HealpixRingProfileWorkspace, 
                                  θ_center, ϕ_center, radius)
            empty!(pixel_list)
            ring_min, ring_max = XGPaint.get_relevant_rings(workspace.res, θ_center, radius)
            
            for ring in ring_min:ring_max
                range1, range2 = XGPaint.get_ring_disc_ranges(workspace, ring, θ_center, ϕ_center, radius)
                
                first_pixel = Healpix.getringinfo(workspace.res, ring).firstPixIdx
                
                for rel_idx in range1; push!(pixel_list, first_pixel + rel_idx - 1); end
                for rel_idx in range2; push!(pixel_list, first_pixel + rel_idx - 1); end
            end
        end
        
        radius = deg2rad(5.0)
        total_pixels = Healpix.nside2npix(nside)
        test_interval = max(1, total_pixels ÷ 50)  # Test every 2% for speed
        
        failures = 0
        for pixel_idx in 1:test_interval:total_pixels
            center_theta, center_phi = Healpix.pix2angRing(res, pixel_idx)
            
            ref_pixels = Healpix.queryDiscRing(res, center_theta, center_phi, radius)
            our_pixels = Int[]
            query_disc_rings!(our_pixels, workspace, center_theta, center_phi, radius)
            
            if Set(ref_pixels) != Set(our_pixels)
                failures += 1
                if failures ≤ 3  # Only report first few failures
                    @warn "Mismatch at pixel $pixel_idx: θ=$(rad2deg(center_theta))°, φ=$(rad2deg(center_phi))°" *
                          " Expected $(length(ref_pixels)), got $(length(our_pixels)) pixels"
                end
            end
        end
        
        @test failures == 0
    end
    
    @testset "Random Points Validation" begin
        nside = 32
        res = Healpix.Resolution(nside)
        workspace = XGPaint.HealpixRingProfileWorkspace(res)
        
        function query_disc_rings!(pixel_list::Vector{Int}, workspace::XGPaint.HealpixRingProfileWorkspace, 
                                  θ_center, ϕ_center, radius)
            empty!(pixel_list)
            ring_min, ring_max = XGPaint.get_relevant_rings(workspace.res, θ_center, radius)
            
            for ring in ring_min:ring_max
                range1, range2 = XGPaint.get_ring_disc_ranges(workspace, ring, θ_center, ϕ_center, radius)
                
                first_pixel = Healpix.getringinfo(workspace.res, ring).firstPixIdx
                
                for rel_idx in range1; push!(pixel_list, first_pixel + rel_idx - 1); end
                for rel_idx in range2; push!(pixel_list, first_pixel + rel_idx - 1); end
            end
        end
        
        import Random
        Random.seed!(42)  # Reproducible results
        test_radius = deg2rad(3.0)
        failures = 0
        
        for i in 1:1000  # Smaller sample for unit test
            center_theta = acos(1 - 2*rand())  # Uniform on sphere
            center_phi = 2π * rand()
            
            ref_pixels = Healpix.queryDiscRing(res, center_theta, center_phi, test_radius)
            our_pixels = Int[]
            query_disc_rings!(our_pixels, workspace, center_theta, center_phi, test_radius)
            
            if Set(ref_pixels) != Set(our_pixels)
                failures += 1
                if failures ≤ 3
                    @warn "Random test $i failed: θ=$(rad2deg(center_theta))°, φ=$(rad2deg(center_phi))°" *
                          " Expected $(length(ref_pixels)), got $(length(our_pixels)) pixels"
                end
            end
        end
        
        @test failures == 0
    end
    
    @testset "Edge Cases" begin
        nside = 32
        res = Healpix.Resolution(nside)
        workspace = XGPaint.HealpixRingProfileWorkspace(res)
        
        # Test exact north pole
        range1, range2 = XGPaint.get_ring_disc_ranges(workspace, 1, 0.0, 0.0, deg2rad(5.0))
        @test !isempty(range1) || !isempty(range2)
        
        # Test exact south pole  
        last_ring = res.nsideTimesFour - 1
        range1, range2 = XGPaint.get_ring_disc_ranges(workspace, last_ring, Float64(π), 0.0, deg2rad(5.0))
        @test !isempty(range1) || !isempty(range2)
        
        # Test very small radius
        range1, range2 = XGPaint.get_ring_disc_ranges(workspace, 50, π/3, π/4, 1e-6)
        # Should return empty or very small ranges
        @test length(range1) + length(range2) ≤ 5
        
        # Test seam crossing (phi near 0/2π)
        range1, range2 = XGPaint.get_ring_disc_ranges(workspace, 50, π/2, 0.1, deg2rad(10.0))
        # Should handle wraparound correctly
        @test range1 isa UnitRange{Int}
        @test range2 isa UnitRange{Int}
    end
    
end
