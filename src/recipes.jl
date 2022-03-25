using Printf
using RecipesBase

function solution_surface( problem :: Problem)

    (η, v, x, t) = solution(problem)
    x, η, t

end

function solution_fourier(problem)

    (η, v, y) = solution(problem)
    fftη = fft(η)
    fftshift(Mesh(y).k), eps(1.0) .+ abs.(fftshift(fftη))

end

@recipe f(::Type{Problem}, problem::Problem) = first(solution(problem))
   

@recipe function f(problem::Problem; fourier = false)

    if fourier

        title --> "Fourier coefficients (log scale)"
        label --> problem.label
        xlabel --> "frequency"
        ylabel --> "amplitude"
        yscale := :log10

        solution_fourier( problem) 

    else

        x, η, t = solution_surface( problem) 

        title --> @sprintf("surface deformation at t= %7.3f", t)
        label --> problem.label
        xlabel --> "x"
        ylabel --> "η"

        x, η

    end

end

@recipe function f(problems::Vector{Problem}; fourier = false)

    layout := (fourier+1, 1)

    for (i,p) in enumerate(problems)

       @series begin
           x, η, t = solution_surface(p)
           title --> @sprintf("Surface deformation at t= %7.3f", t)
           xlabel --> "x"
           ylabel --> "η"
           label --> p.label
           subplot := 1
           x, η
       end

       if fourier
           @series begin
               title --> "Fourier coefficients (log scale)"
               x, y = solution_fourier(p)
               xlabel --> "frequency"
               ylabel --> "amplitude"
               label --> p.label
               subplot := 2
               yscale := :log10
               x, y
           end
       end

    end

end
