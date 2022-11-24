using CirculationModels
using Documenter

DocMeta.setdocmeta!(CirculationModels, :DocTestSetup, :(using CirculationModels); recursive=true)

makedocs(;
    modules=[CirculationModels],
    authors="TS-CUBED <ts-cubed@t-cubed.org.uk> and contributors",
    repo="https://github.com/TS-CUBED/CirculationModels.jl/blob/{commit}{path}#{line}",
    sitename="CirculationModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://TS-CUBED.github.io/CirculationModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Setting up a simple model" => "setupModel.md"
        ],
        "Examples" => [
            "Bjørdalsbakke - Simple Single-Chamber CV-Circuit" => "Bjordalsbakke.md"
        ],
        "Method Index" => "autodoc.md",
    ],
)

deploydocs(;
    repo="github.com/TS-CUBED/CirculationModels.jl",
    devbranch="main",
)





