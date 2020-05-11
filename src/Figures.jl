export create_animation,plot_solution!,plot_solution
using Plots;
using ProgressMeter
"""
	create_animation( p; name::String, ylims=(a,b), Nframes::Int )

Creates an animation showing the evolution of initial-value problems.

Argument `p` is either an element or a vector, list of elements of type `Problem`.

The animation is saved as `name.gif` if `string` is provided.

`ylims` allows to specifies the y axis limits.
If nothing is provided, then the limits are determined from the initial data.
If anything but a Tuple is provided, the axis limits evolve with the solution.

`Nframes` gives the (maximal) number of frames in the animation.

Return `anim` an animation, which can generate a `gif` through `gif(anim, "my_name.gif", fps=15)`

"""
function create_animation( problems; name=nothing, ylims=nothing, Nframes = 201 )
	label=nothing
	if isa(problems,Problem) # if create_animation is called with a single problem
		problems=[problems];  # A one-element array to allow `for pb in p`
		label=problems[1].model.label; # saves the model's name
		problems[1].model.label=""; 	# replace it with blank so that it does not appear in the plots
	end

	if ylims == nothing
		(η,v) = solution(problems[1];t=0)
		y0 = minimum(η);y1 = maximum(η);
		ylims=(y0-(y1-y0)/10,y1+(y1-y0)/10)
	end
	Ns=problems[1].times.Ns
	L=range(1,stop=Ns,step=max(1,ceil(Int,Ns/Nframes)))
	prog = Progress(length(L);dt=1,desc="Creating animation: ")
    anim = @animate for l in L
		plt = plot(layout=(2,1))
		for pb in problems
			plot_solution!( plt, pb; t=pb.times.ts[l] )
		end
		if isa(ylims,Tuple) ylims!(plt[1,1],ylims) end
        next!(prog)
    end
	if label!=nothing problems[1].model.label=label end # puts back the model's name
	if name != nothing gif(anim, string(name,".gif"), fps=15); end
	return anim
end

"""
	plot_solution!( plt; problems; t,x )

Plots in `plt` the solution of initial-value problems at a given time.

Argument `problems` is either an element or a vector or list of elements of type `Problem`.

Keyword `t` is the time. If not provided, then the last computed time is plotted.

If a vector `x` is provided and if it is possible, the solution is interpolated to the collocation points given by `x`.

"""
function plot_solution!( plt, problems; t=nothing,x=nothing )
	if typeof(problems)==Problem problems=[problems] end
	for p in problems
		if x != nothing
			(ηfull,vfull,xfull) = solution(p)
			fftη = fft(ηfull)
			(η,v,X,t) = solution(p;t=t,x=x)
			flag=true
		else
			(η,v,X,t) = solution(p;t=t,x=x)
			fftη = fft(η)
		end

	    plot!(plt[1,1], X, η;
			  title=string("elevation at t=",t),
		      label=p.model.label)

	    plot!(plt[2,1], fftshift(p.mesh.k),
	        log10.(1e-18.+abs.(fftshift(fftη)));
			  title="Fourier modes",
	    	  label=p.model.label)
	end
	nothing
end

"""
	plot_solution( problems; t,x )

Same as `plot_solution!` but generates and returns the plot.
"""
function plot_solution( problems; t=nothing,x=nothing )
	plt = plot(layout=(2,1))
	plot_solution!( plt, problems; t=t,x=x )
	return plt
end

function norm_problem!( plt, p::Problem, s::Real )
	N=[];
	Λ = sqrt.((p.mesh.k.^2).+1);
	prog = Progress(div(p.times.Ns,10),1)
 	for index in range(1,stop=p.times.Ns)
    	(hr,ur) = mapfro(p.model,p.data.U[index])
		push!(N,norm(ifft(Λ.^s.*fft(hr))))
		next!(prog)
	end

    plot!(plt, p.times.ts, N;
		  title=string("norm H^s avec s=",s),
	          label=p.model.label)

end
