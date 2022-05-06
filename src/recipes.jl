using Printf
using RecipesBase

function indices(compression, x) 

    step = compression ? 8 : 1
    1:step:length(x)

end

function solution_surface( problem :: Problem, compression )

    η, v, x, t = solution(problem)
    i = indices(compression, x)
    x[i], η[i], t

end

function solution_surface( problem :: Problem, x̃ :: AbstractArray, compression )

    η, v, x, t = solution(problem; x = x̃, interpolation = false)
    i = indices(compression, x)
    x[i], η[i], t

end

function solution_fourier( problem, compression )

    η, v, y, t = solution(problem)
    fftη = fft(η)
    x = fftshift(Mesh(y).k)
    y = eps(1.0) .+ abs.(fftshift(fftη))
    i = indices(compression, x)
    x[i], y[i]

end

function solution_velocity(problem, compression )

    η, v, x, t = solution(problem)
	x, v, t

end
   

@recipe function f(problem::Problem; var = :surface, compression = false)

	if var == :fourier

		@series begin 

            title --> "Fourier coefficients (log scale)"
            label --> problem.label
            xlabel --> "frequency"
            ylabel --> "amplitude"
            yscale := :log10

            solution_fourier( problem, compression ) 

		end

	end

	if var == :surface

		@series begin 

            x, η, t = solution_surface( problem, compression ) 

            title --> @sprintf("surface deformation at t= %7.3f", t)
            label --> problem.label
            xlabel --> "x"
            ylabel --> "η"

            x, η

		end

	end

	if var == :velocity

		@series begin 

            x, v, t = solution_velocity( problem, compression ) 

            title --> @sprintf("Tangential velocity at t= %7.3f", t)
            label --> problem.label
            xlabel --> "x"
            ylabel --> "v"

            x, v

		end

	end

end

@recipe function f(problem::Problem, x̃::AbstractArray; var = :surface)

    if var == :surface

		@series begin 

            x, η, t = solution_surface( problem, x̃, false )

            title --> @sprintf("surface deformation at t= %7.3f", t)
            label --> problem.label
            xlabel --> "x"
            ylabel --> "η"

            x, η

		end

	end

end


@recipe function f(problems::Vector{Problem}; compression = false)

	surface = get(plotattributes, :surface, true)
	fourier = get(plotattributes, :fourier, false)
	velocity = get(plotattributes, :quiver, false)

    layout := (surface+velocity+fourier, 1)

    for (i,p) in enumerate(problems)

	   n = 0

       if surface
		   n += 1
           @series begin
               x, η, t = solution_surface(p, compression)
               title --> @sprintf("Surface deformation at t=%7.3f", t)
               xlabel --> "x"
               ylabel --> "η"
               label --> p.label
               subplot := n
               x, η
           end
	   end

	   if velocity
		   n += 1
           @series begin
		       x, v, t = solution_velocity(p, compression)
		       title --> @sprintf("Tangential velocity at t=%7.3f",t)
               xlabel --> "x"
               ylabel --> "v"
               subplot := n
		       x, v
	       end
	   end

       if fourier
		   n += 1
           @series begin
               x, y = solution_fourier(p, compression)
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

end



@userplot PlotDifferences

@recipe function f(e::PlotDifferences)

    pairs, problems = _differences_args( e.args )

    for (i, j) in pairs

	    x1, η1, t = solution_surface(problems[i], false)
	    x2, η2, t = solution_surface(problems[j], false)

        @series begin
            xlabel := "x"
            ylabel := "η"
	        label := "$(problems[i].label) - $(problems[j].label)"
	        title := @sprintf("difference (surface deformation) at t=%7.3f", t)
            x1, η1 .- η2
        end

    end

end

function _differences_args( (problem1, problem2 ) :: Tuple{Problem, Problem})
		
     [(1,2)], (problem1, problem2)
           
end

function _differences_args( (problems,) :: Tuple{Vector{Problem}})
		
	@show pairs = [(i,j) for i in eachindex(problems) for j in 1:i-1]
    
    pairs, problems

end
