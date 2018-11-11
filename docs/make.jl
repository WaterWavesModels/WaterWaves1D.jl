push!(LOAD_PATH,"../src/")

using DeepWaterModels
using Documenter
using Plots 
using Literate

# generate example

EXAMPLE     = joinpath(@__DIR__, "..",  "examples/two_problems.jl")
DOC_OUTPUT  = joinpath(@__DIR__, "src", "examples")
NB_OUTPUT   = joinpath(@__DIR__, "src", "notebooks")
Literate.markdown(EXAMPLE, DOC_OUTPUT)
Literate.notebook(EXAMPLE, NB_OUTPUT, execute=false)

EXAMPLE     = joinpath(@__DIR__, "..",  "examples/animation.jl")
DOC_OUTPUT  = joinpath(@__DIR__, "src", "examples")
NB_OUTPUT   = joinpath(@__DIR__, "src", "notebooks")
Literate.markdown(EXAMPLE, DOC_OUTPUT)
Literate.notebook(EXAMPLE, NB_OUTPUT, execute=false)

makedocs(modules=[DeepWaterModels],
         doctest = false,
         format = :html,
         sitename = "DeepWaterModels.jl",
         pages = ["Documentation" => "index.md",
                  "Code basics"   => "basics.md",
                  "Animation"     => "examples/animation.md",
                  "Example"       => "examples/two_problems.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/WaterWavesModels/DeepWaterModels.jl.git",
 )
