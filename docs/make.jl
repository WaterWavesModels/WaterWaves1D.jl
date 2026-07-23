# How-to guides
push!(LOAD_PATH, "../src/")

ENV["GKSwstype"] = "100"
const is_ci = haskey(ENV, "CI")

using Documenter
using DocumenterCitations
using Plots
using Literate
using WaterWaves1D

examples = ["FullDispersion", 
            "HammackSegur", 
            "DeepWater",
            "ShallowWater",
            ]

for example in examples

    INPUT = joinpath(@__DIR__, "..", "examples", example * ".jl")
    OUTPUT = joinpath(@__DIR__, "src", "generated")
    Literate.markdown(INPUT, OUTPUT)

end

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), style = :authoryear)

makedocs(
    modules = [WaterWaves1D],
    plugins = [bib],
    doctest = false,
    authors = "Vincent Duchene and Pierre Navaro",
    format = Documenter.HTML(prettyurls = is_ci),
    sitename = "WaterWaves1D.jl",
    #linkcheck = true,
    #checkdocs = :all,
    warnonly = is_ci ? false : [:cross_references],
    pages = [
        "Home" => "index.md",
        "Background" => "background.md",
        "Main architecture" => "problems.md",
        "How-to..." => "how-to.md",
        "Examples" => ["example.md", ["generated/" * example * ".md" for example in examples]...],
        "Plot recipes" => "plot_recipes.md",
        "Library" => "library.md",
        "References" => "references.md",
    ]
)


deploydocs(
    repo = "github.com/WaterWavesModels/WaterWaves1D.jl.git"
)
