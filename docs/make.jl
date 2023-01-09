using CirculatorySystemModels
using Documenter

DocMeta.setdocmeta!(CirculatorySystemModels, :DocTestSetup, :(using CirculatorySystemModels); recursive=true)

makedocs(;
    modules=[CirculatorySystemModels],
    authors="TS-CUBED <ts-cubed@t-cubed.org.uk> and contributors",
    repo="https://github.com/TS-CUBED/CirculatorySystemModels.jl/blob/{commit}{path}#{line}",
    sitename="CirculatorySystemModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://TS-CUBED.github.io/CirculatorySystemModels.jl",
        edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "Bjørdalsbakke - Simple Single-Chamber CV-Circuit" => "examples/BjordalsbakkeModel.md"
        ],
        "Method Index" => "autodoc.md",
    ]
)

deploydocs(;
    repo="github.com/TS-CUBED/CirculatorySystemModels.jl",
    devbranch="main"
)
