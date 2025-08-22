ENV["GKSwstype"]="100"

using Test
using WaterWaves1D


include("./test_loadsave.jl")
include("./testtypes.jl")
include("./testmodels.jl")
include("./testsolvers.jl")
include("./testinit.jl")
include("./testtools.jl")

include("./test_examples.jl")
include("./test_notebooks.jl")


