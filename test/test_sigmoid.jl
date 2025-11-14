#!/usr/bin/env julia

# Test script to evaluate and plot signmoid τ (optical depth) profiles

using XGPaint
using Plots
using Printf

# ------------------------------------------------------------
# Test parameters - typical cluster
# ------------------------------------------------------------
M_test = 1e14   # solar masses (unitless; multiplied by M_sun internally)
z_test = 0.3

println("Creating τ-profiles for test cluster:")
println("  Mass: $(M_test) M_sun")
println("  Redshift: $(z_test)")

# ------------------------------------------------------------
# Create both standard and sigmoid τ models
# ------------------------------------------------------------
tau_model = BattagliaTauProfile()
sigmoid_tau_model = SigmoidBattagliaTauProfile()

println("Models created successfully!")

# ------------------------------------------------------------
# Compute R_200 and angular size
# ------------------------------------------------------------
R_200_test = XGPaint.R_Δ(tau_model, M_test * XGPaint.M_sun, z_test, 200)
angular_size_test = XGPaint.object_size(tau_model, R_200_test, z_test)

println("Cluster properties:")
println("  R_200 = $(R_200_test / 1000) kpc")  # approximate conversion
println("  Angular size of R_200 = $(angular_size_test * 180/π * 3600) arcsec")

# ------------------------------------------------------------
# Define radial range: 0.01–10 R200 (in angular units)
# ------------------------------------------------------------
r_over_R200 = 10 .^ range(log10(0.01), log10(10), length=100)
radii = r_over_R200 .* angular_size_test   # angular radii in radians

# ------------------------------------------------------------
# Evaluate τ(r)
# ------------------------------------------------------------
println("Evaluating standard Battaglia τ profile...")
tau_standard = [tau_model(r, M_test, z_test) for r in radii]

println("Evaluating sigmoid-suppressed τ profile...")
tau_sigmoid = [sigmoid_tau_model(r, M_test, z_test) for r in radii]

println("Profile evaluation completed!")

# ------------------------------------------------------------
# Plot with Plots.jl (GR backend)
# ------------------------------------------------------------
gr()
plot(r_over_R200, tau_standard;
     lw=2, xscale=:log10, yscale=:log10,
     label="Standard Battaglia τ", color=:blue)
plot!(r_over_R200, tau_sigmoid;
      lw=2, label="Sigmoid-Suppressed τ", color=:red,
      xlabel="r / R₂₀₀", ylabel="Optical Depth τ(r)",
      title="Optical Depth Profiles: Standard vs Sigmoid-Suppressed\nM = 10¹⁴ M☉, z = 0.3")

# Add dashed lines for the sigmoid region
vline!([4.0, 6.0], line=(:dash, :gray), label=["r=4R₂₀₀" "r=6R₂₀₀"])

# Save the plot
savefig("tau_profiles_comparison.png")
println("Plot saved as 'tau_profiles_comparison.png'")

# ------------------------------------------------------------
# Print summary table
# ------------------------------------------------------------
println("\nProfile comparison at key radii:")
println("r/R₂₀₀    Standard τ       Sigmoid τ        Suppression")
println("--------------------------------------------------------")
for r_ratio in [0.1, 1.0, 4.0, 5.0, 6.0, 8.0]
    idx = argmin(abs.(r_over_R200 .- r_ratio))
    t_std = tau_standard[idx]
    t_sig = tau_sigmoid[idx]
    suppression = t_sig / t_std
    @printf("%.1f       %.2e      %.2e      %.3f\n", r_ratio, t_std, t_sig, suppression)
end

# ------------------------------------------------------------
# Check sigmoid suppression function directly
# ------------------------------------------------------------
println("\nSigmoid suppression function values:")
for r_ratio in [3.0, 4.0, 5.0, 6.0, 7.0]
    supp = sigmoid_suppression(r_ratio)
    @printf("sigmoid_suppression(%.1f) = %.4f\n", r_ratio, supp)
end
