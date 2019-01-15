using JLD

## Ne gere pas les param ("NamedTuple"). D'où la traduction en dictionnaire, et vice-versa.
## Pour une raison etrange, ne peut pas enregistrer Px (plan_ifft)
function save_problem(problems::Array{Any,1},name::String)
      problemsave=[]
      for i in range(1,size(problems)[1])
      	push!(problemsave,ProblemSave(problems[i]))
      end

      save(string(name,".jld"),"problemsave",problemsave)
end

function save_problem(problem::Problem,name::String)
      save(string(name,".jld"),"problemsave",ProblemSave(problem))
end

"""
    load_problems(name::String,param::NamedTuple)
    Example of use
    problems=load_problems("foo",param) #"foo.jld" must have bee generated with save_problems

"""
function load_problem(name::String,param::NamedTuple)
      d=load(string(name,".jld"))
      if typeof(d["problemsave"]) == Array{Any,1}
            problems=[]
            for i in range(1,size(d["problemsave"])[1])
            	push!(problems,Problem(d["problemsave"][i],param))
            end
      elseif typeof(d["problemsave"]) == ProblemSave
            problems = Problem(d["problemsave"],param)
      end
      problems
end

function load_problem(name::String)
      d=load(string(name,".jld"))
      if typeof(d["problemsave"]) == Array{Any,1}
            problems=[]
            for i in range(1,size(d["problemsave"])[1])
            	push!(problems,Problem(d["problemsave"][i]))
            end
      elseif typeof(d["problemsave"]) == ProblemSave
            problems = Problem(d["problemsave"])
      end
      problems
end
