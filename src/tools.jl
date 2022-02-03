export interpolate,solution

"""
    interpolate(mesh,vector;n=2^3)

Interpolate a vector `vector` defined on a uniform collocation grid defined by `mesh`.

Returns `(new_mesh,new_vector)` a new uniform mesh with `n` times as many values, and the vector of values on collocation points.

"""
function interpolate(mesh::Mesh,vector;n=2^3::Int)

    fourier=fft(vector)
    m=Int(mesh.N/2)
    new_fourier = [fourier[1:m] ;zeros((n-1)*2*m) ;fourier[m+1:end]]


    new_mesh = Mesh(mesh.xmin,mesh.xmax,n*mesh.N)
    new_vector=ifft(new_fourier)*n
    if all((vector[i] isa Real) for i in length(vector))
        new_vector=real.(new_vector)
    end
    return new_mesh,new_vector
end

"""
    interpolate(mesh,vector,x;fast)

Interpolate a vector `vector` defined on a uniform collocation grid defined by `mesh`, on collocation points given by `x`.

If the collocation points `x` are regularly spaced and the optional keyword argument `fast` is set to `true` (default is `false`),
then the algorithm is faster and uses less allocations, but is less precise.

Returns the vector of values on collocation points.

"""
function interpolate(mesh::Mesh,vector,x;fast=false)

    fourier=fft(vector)
    k = mesh.k
	x₀ = mesh.xmin
	if fast == false
		#same as new_vector=exp.(1im*(x.-x₀)*k')*fourier/length(k)
		new_vector = similar(x,Complex); z = complex.(zero(k));
		for i in 1:length(x)
			z .= exp.(-1im*(x[i].-x₀)*k)
			new_vector[i] = real(dot(z,fourier))
		end
		new_vector ./= length(k)
	else
        if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
            @error("Collocation points must be equally spaced.")
        end
		new_vector = similar(x,Complex); z = exp.(-1im*(x[1].-x₀)*k);eidk=exp.(-1im*(x[2]-x[1])*k);
		for i in 1:length(x)
			new_vector[i] = real(dot(z,fourier))
			z .*= eidk
		end
		new_vector ./= length(k)
	end

    if all((vector[i] isa Real) for i in length(vector))
        new_vector=real.(new_vector)
    end
    return new_vector
end


"""
    solution(p::Problem;t,x,interpolation)

Gives the solution of a solved initial-value at a given time `t`.

# Arguments
- Argument `p` is of type `Problem`.
- Keyword argument `t` is optional, the last computed time is returned by default.
- Keyword argument `x` is optional, if provided the solution is interpolated to the collocation vector `x`.
- Keyword argument `interpolation` is optional, if an integer is provided the solution is interpolated on as many collocation points (if `true`, then the default value `2^3` is chosen).


# Return values
Provides `(η,v,x,t)` where
- `η` is the surface deformation at collocation points;
- `v` is the tangential velocity (derivative of the trace of the velocity potential) at collocation points;
- `x` is the vector of collocation points;
- `t` the time (first computed time greater or equal to provided `t`).

"""
function solution(p::Problem; t=nothing, x=nothing, interpolation = false)
	if isnothing(t) t = p.times.tfin end
	t=min(max(t,0),p.times.tfin)
	index = indexin(false,p.times.ts.<t)[1]
	t=p.times.ts[index]
	if Symbol(typeof(p.model)) == :WaterWaves
		if !isnothing(x) || interpolation != false
			@warn "cannot interpolate non-regularly spaced mesh."
		end
		(x,η,v) = (p.model.mapfro)(p.data.U[index])

	else
    	(η,v) = (p.model.mapfro)(p.data.U[index])
		mesh=p.mesh
		if !isnothing(x)
			η = interpolate(mesh,η,x)
			v = interpolate(mesh,v,x)
		else
			x=mesh.x
		end
		if interpolation == true
			new_mesh,η = interpolate(mesh,η)
			new_mesh,v = interpolate(mesh,v)
			x = new_mesh.x
		elseif isa(interpolation,Int)
			new_mesh,η = interpolate(mesh,η;n=interpolation)
			new_mesh,v = interpolate(mesh,v;n=interpolation)
			x = new_mesh.x
		end
	end
	return η,v,x,t
end
