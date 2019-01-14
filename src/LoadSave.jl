using JLD

## Ne gere pas les param ("NamedTuple"). Chercher comment les traduire en dictionnaire, et vice-versa
## Pour une raison etrange, il semble ne pas pouvoir enregistrer plus que 6 problemes en meme temps...
function save_problem(problems::Array{Problem,1},name::String)
      problemsave=[]
      for i in range(1,size(problems)[1])
      	push!(problemsave,ProblemSave(problems[i]))
      end

      save(string("../../Datas/",name,".jld"),"problemsave",problemsave)
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
      elseif typeof(d["problemsave"]) == Problem
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
      elseif typeof(d["problemsave"]) == Problem
            problems = Problem(d["problemsave"])
      end
      problems
end
