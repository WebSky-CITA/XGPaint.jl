using Documenter, XGPaint

makedocs(;
    modules=[XGPaint],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/xzackli/XGPaint.jl/blob/{commit}{path}#L{line}",
    sitename="XGPaint.jl",
    authors="Zack Li",
    assets=String[],
)

deploydocs(;
    repo="github.com/xzackli/XGPaint.jl",
)
