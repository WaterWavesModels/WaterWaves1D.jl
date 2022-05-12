using Printf
using RecipesBase

function indices(compression, x) 

    step = compression ? 8 : 1
    1:step:length(x)

end

function solution_surface( problem :: Problem, time, compression )

    η, v, x, t = solution(problem, t = time)
    i = indices(compression, x)
    x[i], η[i], t

end

function solution_surface( problem :: Problem, x̃ :: AbstractArray, time, compression )

    η, v, x, t = solution(problem; x = x̃, t = time, interpolation = false)
    i = indices(compression, x)
    x[i], η[i], t

end

function solution_fourier( problem :: Problem, time, compression )

    η, v, y, t = solution(problem, t = time)
    fftη = fft(η)
    x = fftshift(Mesh(y).k)
    y = eps(1.0) .+ abs.(fftshift(fftη))
    i = indices(compression, x)
    x[i], y[i], t

end

function solution_velocity(problem :: Problem, time, compression )

    η, v, x, t = solution(problem, t = time)
    i = indices(compression, x)
    x[i], v[i], t


end

function solution_velocity( problem :: Problem, x̃ :: AbstractArray, time, compression )

    η, v, x, t = solution(problem; x = x̃, t = time, interpolation = false)
    i = indices(compression, x)
    x[i], v[i], t

end
   

@recipe function f(problem::Problem; time = nothing, var = :surface, compression = false)

	titlefontsize --> 10

	if var isa Vector{Symbol}
		variables = var
		layout := (length(var), 1)
	else
        layout := (1, 1)
		variables = [var]
	end

	for (n, variable) in enumerate(variables)

	    if variable == :fourier

	    	@series begin 

                x, y, t = solution_fourier( problem, time, compression ) 

                title --> "Fourier coefficients (log scale)"
                label --> problem.label
                xlabel --> "wavenumber"
                ylabel --> "amplitude"
                yscale := :log10
	    		subplot := n

                x, y

	    	end

	    end

	    if variable == :surface

	    	@series begin 

                x, η, t = solution_surface( problem, time, compression ) 

                title --> @sprintf("surface deformation at t= %7.3f", t)
                label --> problem.label
                xlabel --> "x"
                ylabel --> "η"
	    		subplot := n

                x, η

	    	end

	    end

	    if variable == :velocity

	    	@series begin 

                x, v, t = solution_velocity( problem, time, compression ) 

                title --> @sprintf("Velocity at t= %7.3f", t)
                label --> problem.label
                xlabel --> "x"
                ylabel --> "v"
	    		subplot := n

                x, v

	    	end

	    end

	end

end

@recipe function f(problem::Problem, x̃::AbstractArray; t = nothing, var = :surface)

    if var == :surface

		@series begin 

            x, η, t = solution_surface( problem, x̃, t, false )

            title --> @sprintf("surface deformation at t= %7.3f", t)
            label --> problem.label
            xlabel --> "x"
            ylabel --> "η"

            x, η

		end

	end

end

@recipe function f(problems::Vector{Problem}; var = :surface, time = nothing, compression = false)

	if var isa Vector{Symbol}
		variables = var
		layout := (length(var), 1)
	else
        layout := (1, 1)
		variables = [var]
	end

	layout := (length(variables), 1)
	titlefontsize --> 10

    for (i,p) in enumerate(problems)

	   n = 0

       if :surface in variables
		   n += 1
           @series begin
               x, η, t = solution_surface(p, time, compression)
               title --> @sprintf("Surface deformation at t=%7.3f", t)
               xlabel --> "x"
               ylabel --> "η"
               label --> p.label
               subplot := n
               x, η
           end
	   end

	   if :velocity in variables
		   n += 1
           @series begin
		       x, v, t = solution_velocity(p, time, compression)
		       title --> @sprintf("Velocity at t=%7.3f",t)
               xlabel --> "x"
               ylabel --> "v"
               label --> p.label
               subplot := n
		       x, v
	       end
	   end

       if :fourier in variables
		   n += 1
           @series begin
               x, y, t = solution_fourier(p, time, compression)
               title --> "Fourier coefficients (log scale)"
               xlabel --> "frequency"
               ylabel --> "amplitude"
               label --> p.label
               subplot := n
               yscale := :log10
               x, y
           end
       end

    end

	if :difference in variables || :difference_surface in variables

		x1, η1, t = solution_surface(problems[1], time, compression)
	    x2, η2, t = solution_surface(problems[2], time, compression)

		n += 1

        @series begin
            xlabel := "x"
            ylabel := "Δη"
	        label := "$(problems[1].label) - $(problems[2].label)"
	        title := @sprintf("Difference (surface deformation) at t=%7.3f", t)
            subplot := n
            x1, η1 .- η2
        end

	end

    if :difference_velocity in variables

		x1, v1, t = solution_velocity(problems[1], time, compression)
	    x2, v2, t = solution_velocity(problems[2], time, compression)

		n += 1

        @series begin
            xlabel := "x"
            ylabel := "Δv"
	        label := "$(problems[1].label) - $(problems[2].label)"
	        title := @sprintf("Difference (velocity) at t=%7.3f", t)
            subplot := n
            x1, v1 .- v2
        end

    end

    if :difference_fourier in variables

        x1, y1, t = solution_fourier(problems[1], time, compression)
        x2, y2, t = solution_fourier(problems[2], time, compression)
    
        n += 1
    
        @series begin
            xlabel := "frequency"
            ylabel := "amplitude"
            label := "$(problems[1].label) - $(problems[2].label)"
            title := "Difference (Fourier coefficients in log scale)"
            subplot := n
            yscale := :log10
            x1, abs.(y1 .- y2)
        end

	end

end

@userplot PlotDifferences

@recipe function f(e::PlotDifferences)

    pairs, problems, time = _differences_args( e.args )

    for (i, j) in pairs

	    x1, η1, t = solution_surface(problems[i], time, false)
	    x2, η2, t = solution_surface(problems[j], time, false)

        @series begin
            xlabel := "x"
            ylabel := "η"
	        label := "$(problems[i].label) - $(problems[j].label)"
	        title := @sprintf("Difference (surface deformation) at t=%7.3f", t)
            x1, η1 .- η2
        end

    end

end

function _differences_args( (problem1, problem2 ) :: Tuple{Problem, Problem})
		
     [(1,2)], (problem1, problem2), nothing
           
end

function _differences_args( (problem1, problem2, time ) :: Tuple{Problem, Problem, Real})
		
     [(1,2)], (problem1, problem2), time
           
end

function _differences_args( (problems,) :: Tuple{Vector{Problem}})
		
	@show pairs = [(i,j) for i in eachindex(problems) for j in 1:i-1]
    
    pairs, problems, nothing

end

function _differences_args( (problems, time) :: Tuple{Vector{Problem}, Real})
		
	@show pairs = [(i,j) for i in eachindex(problems) for j in 1:i-1]
    
    pairs, problems, time

end
