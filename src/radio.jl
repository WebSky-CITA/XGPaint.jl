"""Functions for computing radio models."""

abstract type AbstractRadioModel{T<:Real} <: AbstractForegroundModel end

"""
    Radio_Sehgal2009{T}(model parameters...)

Define CIB model parameters. Defaults are from Viero et al. 2013.

```@example
model = CIBModel{Float32}(a_0=0.4)
```
"""
Base.@kwdef struct Radio_Sehgal2009{T<:Real} <: AbstractRadioModel{T}
    nside::Int64    = 4096
    min_redshift::T = 0.0
    max_redshift::T = 4.5
    min_mass::T     = 1e12
    box_size::T     = 40000

    # these coefficients are shared for type I and II
    a_0::T   = 0.0
    a_1_dist::Distributions.Uniform{T} = Distributions.Uniform(T(-0.12), T(0.07))
    a_2_dist::Distributions.Uniform{T} = Distributions.Uniform(T(-0.34), T(0.99))
    a_3_dist::Distributions.Uniform{T} = Distributions.Uniform(T(-0.75), T(-0.23))

    I_R_int::T = 10^(-2.6)
    I_γ::T = 6.0
    I_N_0::T = 1.0
    I_M_0::T = 4e13
    I_α::T = 0.1
    I_L_b::T = 10^(24.0)
    I_m::T = -1.55
    I_n::T = 0.0
    I_δ::T = 3.0
    I_z_p::T = 0.8

    II_R_int::T = 10^(-2.8)
    II_γ::T = 8.0
    II_N_0::T = 0.015
    II_M_0::T = 3e15
    II_α::T = 0.1
    II_L_b::T = 10^(27.5)
    II_m::T = -1.6
    II_n::T = -0.65
    II_z_p::T = 1.3
    II_σ_l::T = 0.4
    II_σ_h::T = 0.73
    I_L_min::T = 2.5e23
    II_L_min::T = 4e23

    # physical constants
    phys_h::T     = 6.62606957e-27      # erg.s
    phys_c::T     = 3e+10               # cm/s
    phys_k::T     = 1.3806488e-16       # erg/K
    phys_Msun::T  = 2e33                # g
    phys_Mpc::T   = 3.086e24            # cm
end

pdf_norm1(x, μ, σ) = exp( -(x-μ)^2 / (2 * σ^2) )

function FR_I_redshift_evolution(z::T, model::Radio_Sehgal2009{T}) where T
    # plateaus at I_z_p
    return (T(1)+min(z, model.I_z_p))^model.I_δ
end

function FR_II_redshift_evolution(z::T, model::Radio_Sehgal2009{T}) where T
    norm = pdf_norm1(T(0), model.II_z_p, model.II_σ_l)
    if z < model.II_z_p
        return pdf_norm1(z, model.II_z_p, model.II_σ_l) / norm
    else
        return pdf_norm1(z, model.II_z_p, model.II_σ_h) / norm
    end
end


export Radio_Sehgal2009
