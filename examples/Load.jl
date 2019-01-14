# # Load previously saved computations
#using DeepWaterModels
include("../src/dependencies.jl")

#----
param = ( ϵ  = 1/2,
        	N  = 2^15,
            L  = 10,
            T  = 5,
            dt = 0.001,
			nr = 100
			)

problemsX1=load_problems("HighFreq1",param)
problemsX2=load_problems("HighFreq2",param)
problemsX3=load_problems("HighFreq3",param)

problemsX4=load_problems("HighFreq4",param)
problemsX5=load_problems("HighFreq5",param)
problemsX6=load_problems("HighFreq6",param)

problems = vcat(problemsX1,problemsX2,problemsX3);
