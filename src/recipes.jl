using Printf
using RecipesBase

function indices(compression :: Bool, x) 

    step = compression ? 8 : 1
    1:step:length(x)

end

function indices(compression :: Int, x) 

    step = compression
    1:step:length(x)

end

function solution_time(problem :: Problem, T)

    η, v, x, time = solution(problem; T = T)
    time

end
function solution_surface( problem :: Problem, T, x̃, interpolation, compression )

    η, v, x = solution(problem; x = x̃, T = T, interpolation = interpolation)
    i = indices(compression, x)
    x[i], η[i]

end

function solution_velocity( problem :: Problem, T, x̃, interpolation, compression )

    η, v, x = solution(problem; x = x̃, T = T, interpolation = interpolation)
    i = indices(compression, x)
    x[i], v[i]

end

function solution_fourier( problem :: Problem, T, x̃, interpolation, compression )

    η, v, x = solution(problem; x = x̃, T = T, interpolation = interpolation)
    fftη = fft(η)
    k = fftshift(Mesh(x).k)
    y = abs.(fftshift(fftη))
    i = indices(compression, x)
    k[i][y[i].!=0], y[i][y[i].!=0]

end

function difference( problems :: Vector{Problem}, T, x̃, interpolation )
    η1, v1, x1 = solution(problems[1]; x = x̃, T = T, interpolation = interpolation)
    η2, v2, x2 = solution(problems[2]; x = x̃, T = T, interpolation = interpolation)

    if !(x1 ≈ x2)

        if Symbol(typeof(problems[2].model)) !== :WaterWaves
            η2, v2, x2, t2 = solution(problems[2]; x = x1, T = T, interpolation = interpolation)

        elseif Symbol(typeof(problems[1].model)) !== :WaterWaves
            η1, v1, x1, t1 = solution(problems[1]; x = x2, T = T, interpolation = interpolation)

        else
            @warn("The difference is computed on different collocation points.")
            η2, v2, x2, t2 = solution(problems[2]; x = x1, T = T, interpolation = interpolation)
        end

    end
    η1, v1, η2, v2, x1

end

function difference_surface( problems :: Vector{Problem}, T, x̃, interpolation, compression )

    η1, v1, η2, v2, x = difference( problems :: Vector{Problem}, T, x̃, interpolation )
    i = indices(compression, x)
    x[i], η1[i].-η2[i]

end

function difference_velocity( problems :: Vector{Problem}, T, x̃, interpolation, compression )

    η1, v1, η2, v2, x = difference( problems :: Vector{Problem}, T, x̃, interpolation )
    i = indices(compression, x)
    x[i], v1[i].-v2[i]

end

function difference_fourier( problems :: Vector{Problem}, T, x̃, interpolation, compression )

    η1, v1, η2, v2, x = difference( problems :: Vector{Problem}, T, x̃, interpolation )
    if (x[2:end].-x[2]≈x[1:end-1].-x[1])
        k = fftshift(Mesh(x).k)
        
    else
        @warn("The FFT is (wrongly) computed on a non-regularly-spaced mesh.")
        k = fftshift(Mesh((xmin=x[1], xmax=x[end]+(x[end]-x[end-1]), N=length(x)) ).k)

    end
    y1 = fftshift(fft(η1))
    y2 = fftshift(fft(η2))
    y= abs.(y1.-y2)

    i = indices(compression, x)
    k[i][y[i].!=0], y[i][y[i].!=0]

end
   

@recipe function f(problem::Problem; 
                        var = :surface, 
                        T = nothing, 
                        x = nothing, 
                        interpolation = false, 
                        compression = false)




	if var isa Vector{Symbol}
		variables = var
	else
		variables = [var]
	end

    layout := (length(variables), 1)
	titlefontsize --> 10


    time = solution_time( problem, T )


	for (n, variable) in enumerate(variables)

        if n == 1 
            string_title = @sprintf(" at t= %7.3f", time)
        else
            string_title = ""
        end

	    if variable == :surface

	    	@series begin 

                title --> string("Surface deformation", string_title)
                label --> problem.label
                xguide --> "x"
                yguide --> "η"
	    		subplot := n

                solution_surface( problem, T, x, interpolation, compression ) 
	    	end

	    end

	    if variable == :velocity

	    	@series begin 


                title --> string("Velocity",string_title)
                label --> problem.label
                xguide --> "x"
                yguide --> "v"
	    		subplot := n

                solution_velocity( problem, T, x, interpolation, compression ) 

	    	end

	    end

        if variable == :fourier || variable == :Fourier

	    	@series begin 


                title --> string("Fourier coefficients (log scale)",string_title)
                label --> problem.label
                xguide --> "wavenumber"
                yguide --> "amplitude"
                yscale --> :log10
	    		subplot := n

                solution_fourier( problem, T, x, interpolation, compression ) 

	    	end

	    end

	end

end

@recipe function f(problems::Vector; 
                    var = :surface, 
                    T = nothing, 
                    x = nothing, 
                    interpolation = false, 
                    compression = false)
    
	if var isa Vector{Symbol}
		variables = var
	else
		variables = [var]
	end

	layout := (length(variables), 1)
	titlefontsize --> 10

    time = solution_time( problems[1], T )

    for (n, variable) in enumerate(variables)

        if n == 1 
            string_title = @sprintf(" at t= %7.3f", time)
        else
            string_title = ""
        end

        if variable == :surface

            for problem in problems

                @series begin 
        
                    title --> string("Surface deformation", string_title)
                    label --> problem.label
                    xguide --> "x"
                    yguide --> "η"
                    subplot := n
    
                    solution_surface( problem, T, x, interpolation, compression ) 
    
                end
    
            end

        end
    
        if variable == :velocity

            for problem in problems
    
                @series begin 
    
    
                    title --> string("Velocity", string_title)
                    label --> problem.label
                    xguide --> "x"
                    yguide --> "v"
                    subplot := n
    
                    solution_velocity( problem, T, x, interpolation, compression ) 
    
                end
    
            end

        end
    
        if variable == :fourier || variable == :Fourier

            for problem in problems
    
                @series begin 
    
    
                    title --> string("Fourier coefficients (log scale)", string_title)
                    label --> problem.label
                    xguide --> "wavenumber"
                    yguide --> "amplitude"
                    yscale --> :log10
                    subplot := n
    
                    solution_fourier( problem, T, x, interpolation, compression ) 

    
                end
    
            end
        end

	    if variable == :difference  || variable == :difference_surface


            @series begin
                xguide --> "x"
                yguide --> "Δη"
	            label --> "$(problems[1].label) - $(problems[2].label)"
	            title --> string("Difference (surface deformation)", string_title)
                subplot := n
                
                difference_surface( [problems[1],problems[2]], T, x, interpolation, compression ) 

            end

	    end

        if variable == :difference_velocity



            @series begin
                xguide --> "x"
                yguide --> "Δv"
                label --> "$(problems[1].label) - $(problems[2].label)"
                title --> string("Difference (velocity)", string_title)
                subplot := n
                
                difference_velocity( [problems[1],problems[2]], T, x, interpolation, compression ) 

            end

        end

        if variable == :difference_fourier || variable == :difference_Fourier


                
            @series begin
                xguide --> "frequency"
                yguide --> "amplitude"
                label --> "$(problems[1].label) - $(problems[2].label)"
                title --> string("Difference (Fourier coefficients in log scale)", string_title)
                yscale --> :log10
                subplot := n
                
                difference_fourier( [problems[1],problems[2]], T, x, interpolation, compression ) 

            end

        end

        if variable == :differences || variable == :differences_surface

            pairs = [(i,j) for i in eachindex(problems) for j in 1:i-1]

            for (i, j) in pairs
                @series begin
                    xguide --> "x"
                    yguide --> "Δη"
                    label --> "$(problems[i].label) - $(problems[j].label)"
                    title --> string("Difference (surface deformation)", string_title)
                    subplot := n
                    
                    difference_surface( [problems[i],problems[j]], T, x, interpolation, compression ) 

                end

            end

        end

        if variable == :differences_velocity

            pairs = [(i,j) for i in eachindex(problems) for j in 1:i-1]

            for (i, j) in pairs

                @series begin
                    xguide --> "x"
                    yguide --> "Δv"
                    label --> "$(problems[i].label) - $(problems[j].label)"
                    title --> string("Difference (velocity)", string_title)
                    subplot := n
                    
                    difference_velocity( [problems[i],problems[j]], T, x, interpolation, compression ) 

                end

            end

        end

        if variable == :differences_fourier || variable == :differences_Fourier

            pairs = [(i,j) for i in eachindex(problems) for j in 1:i-1]

            for (i, j) in pairs
                
                @series begin
                    xguide --> "frequency"
                    yguide --> "amplitude"
                    label --> "$(problems[i].label) - $(problems[j].label)"
                    title --> string("Difference (Fourier coefficients in log scale)", string_title)
                    yscale --> :log10
                    subplot := n
                    
                    difference_fourier( [problems[i],problems[j]], T, x, interpolation, compression ) 

                end

            end

        end


    end

end

@recipe function f(pairs::Vector{Tuple{Problem, Problem}};
                    var = :difference, 
                    T = nothing, 
                    x = nothing, 
                    interpolation = false, 
                    compression = false)


    if var isa Vector{Symbol}
		variables = var
	else
		variables = [var]
	end

	layout := (length(variables), 1)
	titlefontsize --> 10

    for problems in pairs

        time = solution_time( problems[1], T )

        for (n, variable) in enumerate(variables)

            if n == 1 
                string_title = @sprintf(" at t= %7.3f", time)
            else
                string_title = ""
            end

            if variable == :difference  || variable == :difference_surface


                @series begin
                    xguide --> "x"
                    yguide --> "Δη"
                    label --> "$(problems[1].label) - $(problems[2].label)"
                    title --> string("Difference (surface deformation)", string_title)
                    subplot := n
                    
                    difference_surface( [problems[1],problems[2]], T, x, interpolation, compression ) 

                end

            end

            if variable == :difference_velocity



                @series begin
                    xguide --> "x"
                    yguide --> "Δv"
                    label --> "$(problems[1].label) - $(problems[2].label)"
                    title --> string("Difference (velocity)", string_title)
                    subplot := n
                    
                    difference_velocity( [problems[1],problems[2]], T, x, interpolation, compression ) 

                end

            end

            if variable == :difference_fourier || variable == :difference_Fourier


                    
                @series begin
                    xguide --> "frequency"
                    yguide --> "amplitude"
                    label --> "$(problems[1].label) - $(problems[2].label)"
                    title --> string("Difference (Fourier coefficients in log scale)", string_title)
                    yscale --> :log10
                    subplot := n
                    
                    difference_fourier( [problems[1],problems[2]], T, x, interpolation, compression ) 

                end

            end

        end

    end

end