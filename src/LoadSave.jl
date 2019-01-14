using JLD

## Ne gere pas les param ("NamedTuple"). Chercher comment les traduire en dictionnaire, et vice-versa
## Pour une raison etrange, il semble ne pas pouvoir enregistrer plus que 6 problemes en meme temps...
function save_problems(problems::Array{Problem,1},name::String)
      problemsave=[]
      for i in range(1,size(problems)[1])
      	push!(problemsave,ProblemSave(problems[i]))
      end

      save(string("../",name,".jld"),"problemsave",problemsave)
end


"""
    load_problems(problemsave::Array{Any,1},param::NamedTuple)
    Example of use
    @load "foo.jld" #"foo.jld" must have bee generated with save_problems
    problems=load_problems(problemsave,param)

"""
function load_problems(problemsave::Array{Any,1},param::NamedTuple)
      problems=[]
      for i in range(1,size(problemsave)[1])
      	push!(problems,Problem(problemsave[i],param))
      end
      problems
end

"""
    load_problems(name::String,param::NamedTuple)
    Example of use
    problems=load_problems("foo",param) #"foo.jld" must have bee generated with save_problems

"""
function load_problems(name::String,param::NamedTuple)
      d=load(string(name,".jld"))
      problems=[]
      for i in range(1,size(d["problemsave"])[1])
      	push!(problems,Problem(d["problemsave"][i],param))
      end
      problems
end
