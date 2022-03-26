using Printf
using RecipesBase

function solution_surface( problem :: Problem)

    η, v, x, t = solution(problem)
    x, η, t

end

function solution_fourier(problem)

    η, v, y, t = solution(problem)
    fftη = fft(η)
    fftshift(Mesh(y).k), eps(1.0) .+ abs.(fftshift(fftη))

end

function solution_velocity(problem)

    η, v, x, t = solution(problem)
	x, v, t

end
   

@recipe function f(problem::Problem; var = :surface)

	if var == :fourier

		@series begin 

            title --> "Fourier coefficients (log scale)"
            label --> problem.label
            xlabel --> "frequency"
            ylabel --> "amplitude"
            yscale := :log10

            solution_fourier( problem) 

		end

	end

	if var == :surface

		@series begin 

            x, η, t = solution_surface( problem) 

            title --> @sprintf("surface deformation at t= %7.3f", t)
            label --> problem.label
            xlabel --> "x"
            ylabel --> "η"

            x, η

		end

	end

	if var == :velocity

		@series begin 

            x, v, t = solution_velocity( problem) 

            title --> @sprintf("Tangential velocity at t= %7.3f", t)
            label --> problem.label
            xlabel --> "x"
            ylabel --> "v"

            x, v

		end

	end

end

@recipe function f(problems::Vector{Problem})

	surface = get(plotattributes, :surface, true)
	fourier = get(plotattributes, :fourier, false)
	velocity = get(plotattributes, :velocity, false)

	println(velocity)

	delete!(plotattributes, :surface)
	delete!(plotattributes, :velocitx)
	delete!(plotattributes, :fourier)

    layout := (surface+velocity+fourier, 1)

    for (i,p) in enumerate(problems)

	   n = 0

       if surface
		   n += 1
           @series begin
               x, η, t = solution_surface(p)
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
		       x, v, t = solution_velocity(p)
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
               x, y = solution_fourier(p)
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
