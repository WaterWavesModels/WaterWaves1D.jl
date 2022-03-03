# How-to guides
push!(LOAD_PATH,"../src/")

ENV["GKSwstype"]="100"

using Documenter
using Plots
using Literate
using WaterWaves1D

# generate example

DOC_OUTPUT  = joinpath(@__DIR__, "src", "examples")
NB_OUTPUT   = joinpath(@__DIR__, "src", "notebooks")

 #Literate.markdown(joinpath(@__DIR__, "..",  "examples","ManyModelsforWaterWaves.jl"), DOC_OUTPUT, documenter=true)
 #Literate.notebook(joinpath(@__DIR__, "..",  "notebooks","Example.ipynb"), NB_OUTPUT, execute=false)
 #Literate.markdown(joinpath(@__DIR__, "..",  "notebooks","FullDispersion.ipynb"), DOC_OUTPUT, documenter=true)
 #Literate.notebook(joinpath(@__DIR__, "..",  "notebooks","FullDispersion.ipynb"), NB_OUTPUT, execute=false)

makedocs(modules=[WaterWaves1D],
         doctest = false,
         authors = "Vincent Duchene",
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
         sitename = "WaterWaves1D.jl",
         #linkcheck = true,
         #checkdocs = :all,
         pages = [
            "Home" => "home.md",
            "Background" => "background.md",
            "Problems" => "problems.md",
            "How-to..." => "how-to.md",
            "Example" => "example.md",
            "Index" => "index.md"
                   ])

#                  "Animation"     => "examples/animation.md",
#                 "Example"       => "examples/two_problems.md",

deploydocs(
    repo   = "github.com/WaterWavesModels/WaterWaves1D.jl.git"
 )
