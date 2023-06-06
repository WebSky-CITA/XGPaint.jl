using XGPaint
using Documenter
using Healpix
using Cosmology

makedocs(;
    modules=[XGPaint],
    authors="Zack Li",
    repo="https://github.com/WebSky-CITA/XGPaint.jl/blob/{commit}{path}#L{line}",
    sitename="XGPaint.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/WebSky-CITA/XGPaint.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "CIB" => "cib.md",
        "API" => "api.md",
        "Developer Notes" => "developer_notes.md"
    ],
)

deploydocs(;
    repo="github.com/xzackli/XGPaint.jl",
    devbranch = "main"
)
