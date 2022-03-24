using Printf
using RecipesBase

struct WWSolution

    x :: AbstractVector
    y :: AbstractVector
    label :: String
    title :: String
    yscale :: Symbol 

end

function solution_surface( problem :: Problem)

    (η, v, x, t) = solution(problem)
    title = @sprintf("surface deformation at t= %7.3f", t), 
    yscale = :none
    label = problem.label

    WWSolution( x, η, label, title, yscale )

end

function solution_fourier(problem)

    (η, v, y) = solution(problem)
    fftη = fft(η)
    x = fftshift(Mesh(y).k)
    y = eps(1.0) .+ abs.(fftshift(fftη))
    title = "Fourier coefficients (log scale)"
    yscale = :log10
    label = problem.label

    WWSolution( x, y, label, title, yscale )

end


@recipe function f(problem::Problem; fourier = false)

    if fourier
        solution = solution_surface( problem) 
    else
        solution = solution_fourier( problem) 
    end

    x := solution.x
    y := solution.y
    title := solution.title
    label := solution.label
    yscale := solution.yscale
    ()

end
