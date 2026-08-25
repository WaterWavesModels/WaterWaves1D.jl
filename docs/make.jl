# How-to guides
push!(LOAD_PATH, "../src/")

ENV["GKSwstype"] = "100"
const is_ci = haskey(ENV, "CI")

using Documenter
using DocumenterCitations
using DocumenterCodeBlocks
using DocStringExtensions
using Plots
using Literate
using WaterWaves1D

examples = [
    "QuickStart",
    "FullDispersion",
    "HammackSegur",
    "DeepWater",
    "ShallowWater",
]

binder_root_url = "https://plmbinder.math.cnrs.fr/binder/v2/gh/https%3A%2F%2Fplmlab.math.cnrs.fr%2Fnavaro%2FWaterWaves1D.jl/gh-pages?filepath=dev"

for example in examples

    EXAMPLE = joinpath(@__DIR__, "..", "examples", example * ".jl")
    OUTPUT = joinpath(@__DIR__, "src", "generated")
    Literate.markdown(EXAMPLE, OUTPUT; config = Dict("binder_root_url" => binder_root_url))
    Literate.notebook(EXAMPLE, OUTPUT, execute = false)

end

cp(joinpath(@__DIR__, "..", "examples", "Project.toml"), joinpath(@__DIR__, "src", "generated", "Project.toml"); force = true)

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), style = :authoryear)

makedocs(
    modules = [WaterWaves1D],
    plugins = [bib, CodeBlocks()],
    doctest = false,
    authors = "Vincent Duchene and Pierre Navaro",
    format = Documenter.HTML(prettyurls = is_ci),
    sitename = "WaterWaves1D.jl",
    warnonly = is_ci ? false : [:cross_references],
    pages = [
        "Home" => "index.md",
        "Background" => "background.md",
        "Main architecture" => "problems.md",
        "How-to..." => "how-to.md",
        "Examples" => ["generated/" * example * ".md" for example in examples],
        "Plot recipes" => "plot_recipes.md",
        "Library" => "library.md",
        "References" => "references.md",
    ]
)


deploydocs(
    repo = "github.com/WaterWavesModels/WaterWaves1D.jl.git",
    push_preview = true
)
