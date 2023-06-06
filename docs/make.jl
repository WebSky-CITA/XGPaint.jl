using XGPaint
using Documenter
using Healpix
using Cosmology

ENV["PLOTS_DEFAULT_BACKEND"] = "GR"
ENV["GKSwstype"]="nul"
using Plots
using Plots.PlotMeasures: mm

default(
    fontfamily = "Computer Modern", linewidth=1.5,
    titlefontsize=(16+2), guidefontsize=(11+2), 
    tickfontsize=(8+2), legendfontsize=(8+2),
    left_margin=5mm, right_margin=5mm, top_margin = 5mm, bottom_margin = 5mm)

    Plots.default(fontfamily="Computer Modern", fmt=:svg)
makedocs(;
    modules=[XGPaint],
    authors="Zack Li",
    repo="https://github.com/WebSky-CITA/XGPaint.jl/blob/{commit}{path}#L{line}",
    sitename="XGPaint.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/WebSky-CITA/XGPaint.jl",
        assets=["assets/custom.css"]
    ),
    pages=[
        "Home" => "index.md",
        "CIB" => "cib.md",
        "SZ" => "sz.md",
        "API" => "api.md",
        "Developer Notes" => "developer_notes.md"
    ],
)

deploydocs(;
    repo="github.com/xzackli/XGPaint.jl",
    devbranch = "main"
)
