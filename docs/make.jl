push!(LOAD_PATH,"../src/")

using DeepWaterModels
using Documenter
using Plots # to not capture precompilation output

makedocs(modules=[DeepWaterModels],
         doctest = false,
         format = :html,
         sitename = "DeepWaterModels.jl",
         pages = ["Functions" => "index.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/WaterWaveModels/DeepWaterModels.jl.git",
    julia  = "1.0",
    osname = "osx"
 )
