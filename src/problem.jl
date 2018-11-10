export Problem

struct Problem

    model   :: AbstractModel
    initial :: InitialData
    param   :: Parameters
    solver  :: TimeSolver
    data    :: Vector{Tuple{Vector{Complex{Float64}},
			    Vector{Complex{Float64}}}}

    function Problem(model   :: AbstractModel,
         	     initial :: InitialData,
         	     param   :: Parameters,
         	     solver  :: TimeSolver)

         data = [] 

         new(model,initial,param,solver,data)

    end
end

export solve!

function solve!(problem :: Problem, times :: Times) 

    N  = problem.model.mesh.N
    h  = zeros(ComplexF64, N)
    u  = zeros(ComplexF64, N)

    problem.initial(h, u)

    h .= problem.model.Pi .* fft(h)
    u .= problem.model.Pi .* fft(u)
                
    prog = Progress(times.Nt,1) 
    
    push!(problem.data,(h,u))

    for l in range(1,times.Nt-1)
        
        dt = times.t[l+1]-times.t[l]
        
        step!(problem.solver, problem.model, h, u, dt)
    
        push!(problem.data,(h,u))   

        next!(prog)

    end
            
end

