push!(LOAD_PATH,"../src/")

using DeepWaterModels
using Documenter
using Plots # to not capture precompilation output

makedocs(modules=[DeepWaterModels],
         doctest = false,
         format = :html,
         sitename = "DeepWaterModels.jl",
         pages = ["Documentation" => "index.md",
                  "Code basics" => "basics.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/WaterWavesModels/DeepWaterModels.jl.git",
 )
