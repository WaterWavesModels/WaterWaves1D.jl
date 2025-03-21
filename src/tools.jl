export interpolate,solution,mass,momentum,energy,massdiff,momentumdiff,energydiff

"""
    interpolate(mesh,vector;n=2^3)

Interpolate a vector `vector` of values on a uniform collocation grid defined by `mesh`.

Return `(new_mesh,new_vector)` a new uniform mesh with `n` times as many values, and the vector of values at these collocation points.

"""
function interpolate(mesh::Mesh,vector;n=2^3::Int)

    fourier=fft(vector)
    m=Int(mesh.N/2)
    new_fourier = [fourier[1:m] ;zeros((n-1)*2*m) ;fourier[m+1:end]]


    new_mesh = Mesh((xmin=mesh.xmin,xmax=mesh.xmax,N=n*mesh.N))
    new_vector=ifft(new_fourier)*n
    if all((vector[i] isa Real) for i in length(vector))
        new_vector=real.(new_vector)
    end
    return new_mesh,new_vector
end

"""
    interpolate(mesh,vector,x;fast)

Interpolate a vector `vector` of values on a uniform collocation grid defined by `mesh`, on collocation points given by `x`.

If the collocation points `x` are regularly spaced and the optional keyword argument `fast` is set to `true` (default is `false`),
then the algorithm is faster and uses less allocations, but is less precise.

Return the vector of values on collocation points.

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
    solution(pb::Problem;T,x,interpolation)

Give the solution of a solved initial-value problem at a given time `T`.

# Arguments
- Argument `pb` is of type `Problem`.
- Keyword argument `T` is optional, the last computed time is returned by default.
- Keyword argument `x` is optional, if provided the solution is interpolated to the collocation vector `x`.
- Keyword argument `interpolation` is optional, if an integer is provided the solution is interpolated on as many collocation points (if `true`, then the default value `2^3` is chosen).
- Keyword argument `raw` is optional, if set to `true` then `(U,t)` with `U` the raw data and `t` the time is returned (default is `false`).


# Return values
Return `(η,v,x,t)` where
- `η` is the surface deformation at collocation points;
- `v` is the tangential velocity (derivative of the trace of the velocity potential) at collocation points;
- `x` is the vector of collocation points;
- `t` the time (first computed time greater or equal to provided `T`).

"""
function solution(p::Problem; T=nothing, x=nothing, interpolation = false, raw = false)
	if isnothing(T) T = p.times.ts[end] end
	T=min(max(T,0),p.times.ts[end])
	index = findfirst(p.times.ts.>=T)
	t=p.times.ts[index]
	if raw
		return (p.data.U[index]...,t)
	end

	if Symbol(typeof(p.model)) == :WaterWaves || (isnothing(x) && interpolation == false)
		if !isnothing(x) || interpolation != false
			@warn "Cannot interpolate non-regularly spaced mesh."
		end
		return (p.model.mapfro(p.data.U[index])...,t)

	else
		(η,v,y) = p.model.mapfro(p.data.U[index])
		mesh=Mesh(y)
		if isnothing(x)
			x = y
		else
			η = interpolate(mesh,η,x)
			v = interpolate(mesh,v,x)
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
		return η,v,x,t
	end
end

"""
    mass(pb::Problem;T)

Compute the excess of mass of a solved initial-value problem `pb` at a given time `T`.

Keyword argument `T` is optional, the last computed time is used by default.

"""
function mass(p::Problem; T=nothing)
	η,v,x = solution(p;T=T)
	if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
		@error("The excess of mass cannot be computed because the solution is defined on a non-regularly spaced mesh.")
	else
		return sum(η)*(x[2]-x[1])
	end
end

"""
    momentum(pb::Problem;T)

Compute the horizontal impulse of a solved initial-value problem `pb` at a given time `T`.

Keyword argument `T` is optional, the last computed time is used by default.

"""
function momentum(p::Problem; T=nothing)
	η,v,x = solution(p;T=T)
	if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
		@error("The horizontal impulse cannot be computed because the solution is defined on a non-regularly spaced mesh.")
	else
		return sum(η.*v)*(x[2]-x[1])
	end
end

"""
    energy(pb::Problem;T)

Compute the excess of mass of a solved initial-value problem `pb` at a given time `T`.

Keyword argument `T` is optional, the last computed time is used by default.

"""
function energy(p::Problem; T=nothing)
	η,v,x = solution(p;T=T)
	if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
		@error("The energy cannot be computed because the solution is defined on a non-regularly spaced mesh.")
	end
	mesh=Mesh(x)
	U=p.model.mapto(Init(mesh,η,v));
	∂ₓ=1im*mesh.k;
	p.model.f!(U);fftm=-U[1]./∂ₓ;fftm[1]=0;
	m=real(ifft(fftm));
	@. mesh.dx/2*($sum(η^2) + $sum(v*m))
end

"""
    massdiff(pb::Problem;T,rel)

Compute the difference of excess of mass of a solved initial-value problem `pb` between given time `T` and initial time.

Keyword argument `T` is optional, the last computed time is used by default.

If keyword argument `rel=true` (default is false), then compute the relative difference (with initial value as reference).

"""
function massdiff(p::Problem; T=nothing,rel=false)
	η,v,x = solution(p;T=T)
	η0,v0,x0 = solution(p;T=0)
	if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
		@error("The excess of mass difference cannot be computed because the solution is defined on a non-regularly spaced mesh.")
	else
		if rel==false return sum(η-η0)*(x[2]-x[1]) else return sum(η-η0)/sum(η0) end
	end
end

"""
    momentumdiff(pb::Problem;T,rel)

Compute the difference of horizontal impulse of a solved initial-value problem `pb` between given time `T` and initial time.

Keyword argument `T` is optional, the last computed time is used by default.

If keyword argument `rel=true` (default is false), then compute the relative difference (with initial value as reference).

"""
function momentumdiff(p::Problem; T=nothing,rel=false)
	η,v,x = solution(p;T=T)
	η0,v0,x0 = solution(p;T=0)
	if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
		@error("The horizontal impulse difference cannot be computed because the solution is defined on a non-regularly spaced mesh.")
	else
		if rel==false return sum((η-η0).*v+η0.*(v-v0))*(x[2]-x[1]) else return sum((η-η0).*v+η0.*(v-v0))/sum(η0.*v0) end
	end
end

"""
    energydiff(pb::Problem;T,rel)

Compute the difference of energy of a solved initial-value problem `pb` between given time `T` and initial time.

Keyword argument `T` is optional, the last computed time is used by default.

If keyword argument `rel=true` (default is false), then compute the relative difference (with initial value as reference).

"""
function energydiff(p::Problem; T=nothing, rel=false)
	η,v,x = solution(p;T=T)
	η0,v0,x0 = solution(p;T=0)
	if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
		@error("The energy difference cannot be computed because the solution is defined on a non-regularly spaced mesh.")
	end
	mesh=Mesh(x)
	U=p.model.mapto(Init(mesh,η,v)); 
	U0=p.model.mapto(Init(mesh,η0,v0));
	∂ₓ=1im*mesh.k;
	p.model.f!(U);fftm=-U[1]./∂ₓ;fftm[1]=0;
	p.model.f!(U0);fftm0=-U0[1]./∂ₓ;fftm0[1]=0;
	m=real(ifft(fftm));	m0=real(ifft(fftm0));
	if rel == false
		return @. mesh.dx/2*($sum((η-η0)*η+η0*(η-η0)) + $sum((v-v0)*m+v0*(m-m0)))
	else
		return @. ($sum((η-η0)*η+η0*(η-η0)) + $sum((v-v0)*m+v0*(m-m0)))/($sum(η0^2) + $sum(v0*m0))
	end
end
