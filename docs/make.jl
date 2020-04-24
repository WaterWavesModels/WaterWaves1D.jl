push!(LOAD_PATH,"../src/")

ENV["GKSwstype"]="100"

using ShallowWaterModels
using Documenter
using Plots 
using Literate

# generate example

DOC_OUTPUT  = joinpath(@__DIR__, "src", "examples")
NB_OUTPUT   = joinpath(@__DIR__, "src", "notebooks")

Literate.markdown(joinpath(@__DIR__, "..",  "examples","two_problems.jl"), DOC_OUTPUT, documenter=true)
Literate.notebook(joinpath(@__DIR__, "..",  "examples","two_problems.jl"), NB_OUTPUT, execute=false)
Literate.markdown(joinpath(@__DIR__, "..",  "examples","animation.jl"), DOC_OUTPUT, documenter=true)
Literate.notebook(joinpath(@__DIR__, "..",  "examples","animation.jl"), NB_OUTPUT, execute=false)

makedocs(modules=[ShallowWaterModels],
         doctest = false,
         authors = "Vincent Duchene",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         sitename = "ShallowWaterModels.jl",
         pages = ["Documentation" => "index.md",
                  "Code basics"   => "basics.md",
                  "Animation"     => "examples/animation.md",
                  "Example"       => "examples/two_problems.md",
                  "Contents"      => "contents.md"])

deploydocs(
    repo   = "github.com/WaterWavesModels/ShallowWaterModels.jl.git"
 )
