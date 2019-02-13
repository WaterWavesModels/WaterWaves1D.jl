import JLD

export save, load

## Ne gere pas les param ("NamedTuple"). D'où la traduction en dictionnaire, et vice-versa.
## Pour une raison etrange, ne peut pas enregistrer Px (plan_ifft)
function save(problems :: Vector{Problem}, name::String)

    problemsave = ProblemSave[]
    for p in problems
        push!(problemsave,convert(ProblemSave,p))
    end

    JLD.save(string(name,".jld"),"problemsave",problemsave)

end

function save(problem::Problem,name::String)
    JLD.save(string(name,".jld"),"problemsave",convert(ProblemSave,problem))
end

"""
    load_problems(name::String,param::NamedTuple)
    Example of use
    problems=load_problems("foo",param) #"foo.jld" must have bee generated with save_problems

"""
function load(name::String)

      problems = Problem[]

      for p in JLD.load(string(name,".jld"))
          push!(convert(Problem,p), problems)
      end

      problems

end
