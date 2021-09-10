push!(LOAD_PATH,"../src/")

ENV["GKSwstype"]="100"

using Documenter
using Plots 
using Literate
using WaterWaves1D

# generate example

DOC_OUTPUT  = joinpath(@__DIR__, "src", "examples")
NB_OUTPUT   = joinpath(@__DIR__, "src", "notebooks")

# Literate.markdown(joinpath(@__DIR__, "..",  "examples","two_problems.jl"), DOC_OUTPUT, documenter=true)
# Literate.notebook(joinpath(@__DIR__, "..",  "examples","two_problems.jl"), NB_OUTPUT, execute=false)
# Literate.markdown(joinpath(@__DIR__, "..",  "examples","animation.jl"), DOC_OUTPUT, documenter=true)
# Literate.notebook(joinpath(@__DIR__, "..",  "examples","animation.jl"), NB_OUTPUT, execute=false)

makedocs(modules=[WaterWaves1D],
         doctest = false,
         authors = "Vincent Duchene",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         sitename = "WaterWaves1D.jl",
         pages = [ "Quickstart" => "quickstart.md", 
                   "Documentation" => "index.md",
                   "Code basics"   => "basics.md",
                   "Contents"      => "contents.md"])

#                  "Animation"     => "examples/animation.md",
#                 "Example"       => "examples/two_problems.md",

deploydocs(
    repo   = "github.com/WaterWavesModels/WaterWaves1D.jl.git"
 )
