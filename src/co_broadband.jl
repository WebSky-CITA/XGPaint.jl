"""Functions for computing broadband CO/[C I] models."""

using SparseArrays

abstract type AbstractCOModel{T<:Real} <: AbstractForegroundModel end

"""
    CO_CROWNED{T}(; kwargs...)

Define broadband CO model parameters. All numbers not typed are converted to type T. This model has the following parameters and default values:

* `nside::Int64 = 4096`
* `hod::String = "shang"`

"""
@with_kw struct CO_CROWNED{T<:Real} <: AbstractCOModel{T} @deftype T
    nside::Int64    = 4096
    Jupmax::Int64   = 7
    alpha_IR_CO     = 1.37
    beta_IR_CO      =-1.74
    scatterdex_IR_CO= 0.3
    r21_average     = 0.75
    r21_sigmascatter= 0.11
    r21_rolloff     = 0.86 # 1 cuts off r21 sharply; set < 1 to avoid errors
    R_CO10_CI       = 0.18*77.83 # r=0.18 times cubic frequency scaling
    scatterdex_CO_CI= 0.2
    xRJ_GHz = 0.0176086 # factor used in RJ to thermodynamic deltaT conversion
end

# fit at 2.10833 THz good within 2% for 20-21.5 K, within 0.5% for 21.5-50K
Lnu_to_LIR(Td) = 9.521f8*(53.865 + 3.086*Td - 0.213*Td^2 + 0.00749*Td^3);

# process sources NamedTuple from AbstractCIBModel generate_sources
function process_sources(model::AbstractCOModel{T}, sources,
            cib_model::AbstractCIBModel{T}; verbose=true) where T
    r21_ratio = (1.0f0-model.r21_average)/model.r21_sigmascatter
    r21_ratioA = r21_ratio*model.r21_rolloff
    r21_ratioB = r21_ratio*(2-model.r21_rolloff)
    verbose && println("Allocating CO arrays for $(sources.N_cen) centrals.")
    LIR_cen = Array{Float32,1}(undef,sources.N_cen)
    LcoJ_cen = Array{Float32,2}(undef,sources.N_cen,model.Jupmax)
    nuJ_cen = Array{Float32,2}(undef,(sources.N_cen,model.Jupmax))
    quasiTcoJ_cen = Array{Float32,2}(undef,(sources.N_cen,model.Jupmax))
    nuCI_cen = Array{Float32,1}(undef,sources.N_cen)
    LCI_cen = Array{Float32,1}(undef,sources.N_cen)
    quasiTCI_cen = Array{Float32,1}(undef,sources.N_cen)
    verbose && println("Allocating CO arrays for $(sources.N_sat) centrals.")
    LIR_sat = Array{Float32,1}(undef,sources.N_sat)
    LcoJ_sat = Array{Float32,2}(undef,sources.N_sat,model.Jupmax)
    nuJ_sat = Array{Float32,2}(undef,(sources.N_sat,model.Jupmax))
    quasiTcoJ_sat = Array{Float32,2}(undef,(sources.N_sat,model.Jupmax))
    nuCI_sat = Array{Float32,1}(undef,sources.N_sat)
    LCI_sat = Array{Float32,1}(undef,sources.N_sat)
    quasiTCI_sat = Array{Float32,1}(undef,sources.N_sat)
    Threads.@threads :static for i in 1:sources.N_cen
        # calculate dust temperatures and bolometric LIR for each source
        Td = cib_model.shang_Td * (1.f0 .+sources.redshift_cen[i])^cib_model.shang_alpha;
        LIR_cen[i] = Lnu_to_LIR(Td)*sources.lum_cen[i]*nu2theta(2.10833f12,sources.redshift_cen[i],cib_model)
        rnd = Float32(randn())
        r21 = model.r21_average+model.r21_sigmascatter*(rnd+r21_ratioB-sqrt(rnd^2.0f0-2*r21_ratioA*rnd+r21_ratioB^2))/2.0f0 # cuts off r21 at 1, assumes subthermal
        # generate CO SLED for each source with log-normal LCO(1-0) scatter
        LcoJ_cen[i,1] = LIR_cen[i]^(1/model.alpha_IR_CO)*10^(-model.beta_IR_CO/model.alpha_IR_CO)*exp((randn()-0.5*2.302585*model.scatterdex_IR_CO)*2.302585*model.scatterdex_IR_CO)*4.9e-5
        for J in 2:model.Jupmax
            LcoJ_cen[i,J] = LcoJ_cen[i,1]/(1.0f0+exp((log.(1.0f0/r21-1.0f0)-1.0f0)+Float32(J-1)))*J^3
        end
        for J in 1:model.Jupmax
            nuJ_cen[i,J] = 115.27f0*J/(1.0f0+sources.redshift_cen[i])
            quasiTcoJ_cen[i,J] = LcoJ_cen[i,J]*(1.315415f0/4.0f0/pi/nuJ_cen[i,J]^2.0f0/sources.dist_cen[i]^2.0f0/(1.0f0+sources.redshift_cen[i])^2.0f0) # divide by dnu in GHz
            xRJ = model.xRJ_GHz*nuJ_cen[i,J]
            quasiTcoJ_cen[i,J]/= xRJ^2*exp(xRJ)/(exp(xRJ)-1)^2
        end
        nuCI_cen[i] = 492.16f0/(1.0f0+sources.redshift_cen[i])
        LCI_cen[i] = LcoJ_cen[i,1]*model.R_CO10_CI*exp((randn()-0.5*2.302585*model.scatterdex_CO_CI)*2.302585*model.scatterdex_CO_CI)
        quasiTCI_cen[i] = LCI_cen[i]*(1.315415f0/4.0f0/pi/nuCI_cen[i]^2.0f0/sources.dist_cen[i]^2.0f0/(1.0f0+sources.redshift_cen[i])^2.0f0) # divide by dnu in GHz
        xRJ = model.xRJ_GHz*nuCI_cen[i]
        quasiTCI_cen[i]/= xRJ^2*exp(xRJ)/(exp(xRJ)-1)^2
    end
    Threads.@threads :static for i in 1:sources.N_sat
        Td = cib_model.shang_Td * (1.f0 .+sources.redshift_sat[i])^cib_model.shang_alpha;
        LIR_sat[i] = Lnu_to_LIR(Td)*sources.lum_sat[i]*XGPaint.nu2theta(2.10833f12,sources.redshift_sat[i],cib_model)
        rnd = Float32(randn())
        r21 = model.r21_average+model.r21_sigmascatter*(rnd+r21_ratioB-sqrt(rnd^2.0f0-2*r21_ratioA*rnd+r21_ratioB^2))/2.0f0 # cuts off r21 at 1, assumes subthermal
        # generate CO SLED for each source with log-normal LCO(1-0) scatter
        LcoJ_sat[i,1] = LIR_sat[i]^(1/model.alpha_IR_CO)*10^(-model.beta_IR_CO/model.alpha_IR_CO)*exp((randn()-0.5*2.302585*model.scatterdex_IR_CO)*2.302585*model.scatterdex_IR_CO)*4.9e-5
        LcoJ_sat[i,1] = LIR_sat[i]^(1/1.37)*10^(1.74/1.37)*exp((randn()-0.5*2.302585*0.3)*2.302585*0.3)*4.9e-5
        for J in 2:model.Jupmax
            LcoJ_sat[i,J] = LcoJ_sat[i,1]/(1.0f0+exp((log.(1.0f0/r21-1.0f0)-1.0f0)+Float32(J-1)))*J^3
        end
        for J in 1:model.Jupmax
            nuJ_sat[i,J] = 115.27f0*J/(1.0f0+sources.redshift_sat[i])
            quasiTcoJ_sat[i,J] = LcoJ_sat[i,J]*(1.315415f0/4.0f0/pi/nuJ_sat[i,J]^2.0f0/sources.dist_sat[i]^2.0f0/(1.0f0+sources.redshift_sat[i])^2.0f0) # divide by dnu in GHz
            xRJ = model.xRJ_GHz*nuJ_sat[i,J]
            quasiTcoJ_sat[i,J]/= xRJ^2*exp(xRJ)/(exp(xRJ)-1)^2
        end
        nuCI_sat[i] = 492.16f0/(1.0f0+sources.redshift_sat[i])
        LCI_sat[i] = LcoJ_sat[i,1]*model.R_CO10_CI*exp((randn()-0.5*2.302585*model.scatterdex_CO_CI)*2.302585*model.scatterdex_CO_CI)
        quasiTCI_sat[i] = LCI_sat[i]*(1.315415f0/4.0f0/pi/nuCI_sat[i]^2.0f0/sources.dist_sat[i]^2.0f0/(1.0f0+sources.redshift_sat[i])^2.0f0) # divide by dnu in GHz
        xRJ = model.xRJ_GHz*nuCI_sat[i]
        quasiTCI_sat[i]/= xRJ^2*exp(xRJ)/(exp(xRJ)-1)^2
    end
    return (
        N_cen = sources.N_cen,
        N_sat = sources.N_sat,
        hp_ind_cen = sources.hp_ind_cen,
        hp_ind_sat = sources.hp_ind_sat,
        LIR_cen = LIR_cen,
        LcoJ_cen = LcoJ_cen,
        nuJ_cen = nuJ_cen,
        quasiTcoJ_cen = quasiTcoJ_cen,
        nuCI_cen = nuCI_cen,
        LCI_cen = LCI_cen,
        quasiTCI_cen = quasiTCI_cen,
        LIR_sat = LIR_sat,
        LcoJ_sat = LcoJ_sat,
        nuJ_sat = nuJ_sat,
        quasiTcoJ_sat = quasiTcoJ_sat,
        nuCI_sat = nuCI_sat,
        LCI_sat = LCI_sat,
        quasiTCI_sat = quasiTCI_sat
    )
end

# paint for co_broadband takes bandpass
#     bandpass is a Nx2 array with column 1 giving central frequencies (GHz)
#                                  column 2 giving arbitrary transmission
function paint!(result_map_CO::HealpixMap{T_map, RingOrder},
        result_map_CI::HealpixMap{T_map, RingOrder},
        bandpass::AbstractArray,
        model::AbstractCOModel{T}, sources) where {T_map, T}
    fluxCO_cen = Array{Float32,2}(undef,(sources.N_cen,model.Jupmax))
    fluxCO_sat = Array{Float32,2}(undef,(sources.N_sat,model.Jupmax))
    fluxCI_cen = Array{Float32,1}(undef,sources.N_cen)
    fluxCI_sat = Array{Float32,1}(undef,sources.N_sat)
    bandpass_edges = [bandpass[1,1].-diff(bandpass[1:2,1])/2;bandpass[1:end-1,1]+diff(bandpass[:,1])/2;bandpass[end,1].+diff(bandpass[end-1:end,1])/2]
    bandpass_norm = bandpass[:,2]/sum(bandpass[:,2])
    bandpass_diff = diff(bandpass_edges)
    bandpass_len = length(bandpass_norm)
    mCO = result_map_CO.pixels
    mCO_pixsr = nside2pixarea(result_map_CO.resolution.nside)
    fill!(mCO, zero(T))
    mCI = result_map_CI.pixels
    mCI_pixsr = nside2pixarea(result_map_CI.resolution.nside)
    fill!(mCI, zero(T))
    println("painting CO/CI from ",sources.N_cen," centrals")
    Threads.@threads :static for i in 1:sources.N_cen
        for J in 1:model.Jupmax
            j = searchsortedlast(bandpass_edges,sources.nuJ_cen[i,J])
            if (j>0) && (j<=bandpass_len)
                fluxCO_cen[i,J] = sources.quasiTcoJ_cen[i,J]*bandpass_norm[j]/bandpass_diff[j]/mCO_pixsr;
                mCO[sources.hp_ind_cen[i]]+= fluxCO_cen[i,J];
            end
        end
        j = searchsortedlast(bandpass_edges,sources.nuCI_cen[i])
        if (j>0) && (j<=bandpass_len)
            fluxCI_cen[i] = sources.quasiTCI_cen[i]*bandpass_norm[j]/bandpass_diff[j]/mCO_pixsr;
            mCI[sources.hp_ind_cen[i]]+= fluxCI_cen[i];
        end
    end
    println("painting CO/CI from ",sources.N_sat," satellites")
    Threads.@threads :static for i in 1:sources.N_sat
        for J in 1:model.Jupmax
            j = searchsortedlast(bandpass_edges,sources.nuJ_sat[i,J])
            if (j>0) && (j<=bandpass_len)
                fluxCO_sat[i,J] = sources.quasiTcoJ_sat[i,J]*bandpass_norm[j]/bandpass_diff[j]/mCO_pixsr;
                mCO[sources.hp_ind_sat[i]]+= fluxCO_sat[i,J];
            end
        end
        j = searchsortedlast(bandpass_edges,sources.nuCI_sat[i])
        if (j>0) && (j<=bandpass_len)
            fluxCI_sat[i] = sources.quasiTCI_sat[i]*bandpass_norm[j]/bandpass_diff[j]/mCO_pixsr;
            mCI[sources.hp_ind_sat[i]]+= fluxCI_sat[i];
        end
    end
    return fluxCO_cen, fluxCI_cen, fluxCO_sat, fluxCI_sat
end
