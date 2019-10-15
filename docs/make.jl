using Documenter, XGPaint

makedocs(;
    modules=[XGPaint],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/xzackli/XGPaint/blob/{commit}{path}#L{line}",
    sitename="XGPaint",
    authors="Zack Li",
    assets=String[],
)

deploydocs(;
    repo="github.com/xzackli/XGPaint",
)
