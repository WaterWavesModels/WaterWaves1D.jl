using NFFT
export myifft, interpolate, solution, mass, momentum, energy, mass_diff, momentum_diff, energy_diff

"""
$(TYPEDSIGNATURES)

Inverse Fast Fourier transform applied to the vector of fourier coefficients `fourier` of a `P`-periodic signal, 
evaluated at non-uniformly distributed collocation points `sampling`.

A relative error tolerance can be specified with the keyword argument `reltol` (default value is `reltol=1e-8`).

A translation value `x₀` may be specified (by default, `x₀`, that is the uniform grid defining Fourier coefficients is centered).


Return the vector of values at collocation points `sampling`.
"""
function myifft(fourier , sampling, P ; x₀ = 0, reltol=1e-8)
    nfft(  -mod.((sampling.-x₀)/P,1) .+ 1/2,ifftshift(fourier),reltol=reltol)/length(fourier)
end


"""
$(TYPEDSIGNATURES)

Interpolate a vector `vector` of values on a uniform collocation grid defined by `mesh`.

Return `(new_mesh,new_vector)` a new uniform mesh with `n` times as many values, and the vector of values at these collocation points.

"""
function interpolate(mesh::Mesh, vector; n = 2^3::Int)

    fourier = fft(vector)
    m = Int(mesh.N / 2)
    new_fourier = [fourier[1:m] ;zeros((n - 1) * 2 * m) ;fourier[(m + 1):end]]


    new_mesh = Mesh((xmin = mesh.xmin, xmax = mesh.xmax, N = n * mesh.N))
    new_vector = ifft(new_fourier) * n
    if all((vector[i] isa Real) for i in length(vector))
        new_vector = real.(new_vector)
    end
    return new_mesh, new_vector
end

"""
$(TYPEDSIGNATURES) 

Interpolate a vector `vector` of values on a uniform collocation grid defined by `mesh`, on collocation points given by `x`.

If the optional keyword argument `fast` is set to `true` (default is `false`),
then the algorithm is faster and uses less allocations, but may be less precise 
(a relative error estimation may be provided by the keyword argument `reltol`; default value is `reltol=1e-8`).

Return the vector of values on collocation points.

"""
function interpolate(mesh::Mesh, vector, x; fast = false, reltol = 1e-8)

    fourier = fft(vector)
    k = mesh.k
    x₀ = mesh.xmin
    if fast == false
        #same as new_vector=exp.(1im*(x.-x₀)*k')*fourier/length(k)
        new_vector = similar(x, Complex); z = complex.(zero(k))
        for i in 1:length(x)
            z .= exp.(-1im * (x[i] .- x₀) * k)
            new_vector[i] = real(dot(z, fourier))
        end
        new_vector ./= length(k)
    else
        new_vector = myifft(fourier,x,mesh.xmax-mesh.xmin; x₀=(mesh.xmax+mesh.xmin)/2, reltol=reltol)
    end

    if all((vector[i] isa Real) for i in length(vector))
        new_vector = real.(new_vector)
    end
    return new_vector
end

"""
$(TYPEDSIGNATURES) 

Interpolate a vector `vector` of values on a non-uniform grid `y`, on collocation points given by `x`.

If the grid `y` is uniform, consider other interpolation methods using meshes.

If the optional keyword argument `fast` is set to `true` (default is `false`),
then the algorithm is faster and uses less allocations, but may be less precise 
(a relative error estimation may be provided by the keyword argument `reltol`; default value is `reltol=1e-8`).

Return the vector of values on collocation points.

"""
function interpolate(y, vector, x; P=nothing, fast = false, reltol=1e-8)

    # If not provided by the user, infer (approximate) half-period from the values of the non-uniform grid
    if isnothing(P) 
        P=(y[end]-y[1])*(1+1/length(y))
    end
    mesh=Mesh((L=P/2,N=length(y))) # uniform mesh with same size as y
 
   # We seek Fourier frequencies f such that the associated function has values 'vector' on the points 'y'
    
   if fast == true  # iterative method (works if y is near-uniform)
        Π=abs.(mesh.k).<1/maximum(abs.(y-mesh.x).+eps()); # frequency cut-off
        f0=fft(vector); # initial guess
        f=copy(f0);
        δf=copy(f0);
        
        flag = false; iter=0;
        rel_error = maximum(abs.(vector-myifft(f,y,P;x₀=y[1] +P/2,  reltol=reltol/10)))/maximum(abs.(vector)) 
        new_rel_error = 0
        while iter < 20 && rel_error > reltol
            iter+=1
            δf.=f0-fft(myifft(f0,y,P;x₀=y[1] +P/2,reltol=reltol/10)) # compute next term in the asymptotic expansion
            f.+=Π.*δf # add the new term
            f0.=Π.*δf # prepare for next iteration
            new_rel_error=maximum(abs.(vector-myifft(f,y,P;x₀=y[1] +P/2,reltol=reltol)))/maximum(abs.(vector))
            if new_rel_error > rel_error
                f.-=δf # undo the last step
                flag = true
                break # stop iteration
            else
                rel_error=new_rel_error
            end
        end
        if flag || iter == 20
            @warn("relative error target was not achieved after $iter iterations: $rel_error.
            Consider using optional argument 'fast=false'")
        end

    else # Linear algebra method
        f=exp.(1im*(y.-mesh.xmin.-y[1] .-P/2)*mesh.k')/length(mesh.k) \ vector
        rel_error=maximum(abs.(vector-myifft(f,y,P;x₀=y[1] +P/2,reltol=reltol)))/maximum(abs.(vector))
    end
    @warn("relative error estimation: $rel_error")
    new_vector = myifft(f,x,P;x₀=y[1] +P/2,reltol=reltol)

    if all((vector[i] isa Real) for i in length(vector))
        new_vector = real.(new_vector)
    end

    return new_vector

end


"""
$(TYPEDSIGNATURES)

Give the solution of a solved initial-value problem at a given time `T`.

# Arguments
- Argument `pb` is of type `Problem`.
- Keyword argument `T` is optional, the last computed time is returned by default.
- Keyword argument `x` is optional, if provided the solution is interpolated to the vector of collocation points `x`.
- Keyword arguments `fast` and `reltol` are arguments of the `interpolate` function (by default `fast=false` and `reltol=1e-8`).
- Keyword argument `raw` is optional, if set to `true` then `(U,t)` with `U` the raw data and `t` the time is returned (default is `false`).


# Return values
Return `(η,v,x,t)` where
- `η` is the surface deformation at collocation points;
- `v` is the tangential velocity (derivative of the trace of the velocity potential) at collocation points;
- `x` is the vector of collocation points;
- `t` the time (first computed time greater or equal to provided `T`).

"""
function solution(p::Problem; T = nothing, x = nothing, fast = false, reltol=1e-8, raw = false)
    if isnothing(T)
        T = p.times.ts[end]
    end
    T = min(max(T, 0), p.times.ts[end])
    index = findfirst(p.times.ts .>= T)
    t = p.times.ts[index]
    if raw
        return (p.data.U[index]..., t)
    end
    if isnothing(x)
        return (p.model.mapfro(p.data.U[index])...,t)
    else # interpolation on vector of collocation points `x`.
        (η, v, y) = p.model.mapfro(p.data.U[index])
        if (y[2:end] .- y[2] ≈ y[1:(end - 1)] .- y[1]) # equally spaced collocation points
            mesh = Mesh(y)
            η = interpolate(mesh, η, x)
            v = interpolate(mesh, v, x)
        else
            η = interpolate(y, η, x)
            v = interpolate(y, v, x)
        end
        return η, v, x, t
    end
    
end

"""
$(TYPEDSIGNATURES)

Compute the excess of mass of a solved initial-value problem `pb` at a given time `T`.

Keyword argument `T` is optional, the last computed time is used by default.

"""
function mass(p::Problem; T = nothing)::Float64
    η, v, x = solution(p; T = T)

    @assert (x[2:end] .- x[2] ≈ x[1:(end - 1)] .- x[1]) "The excess of mass cannot be computed because the solution is defined on a non-regularly spaced mesh."

    return sum(η) * (x[2] - x[1])
end

"""
$(TYPEDSIGNATURES)

Compute the horizontal impulse of a solved initial-value problem `pb` at a given time `T`.

Keyword argument `T` is optional, the last computed time is used by default.

"""
function momentum(p::Problem; T = nothing)
    η, v, x = solution(p; T = T)
    if !(x[2:end] .- x[2] ≈ x[1:(end - 1)] .- x[1])
        @error("The horizontal impulse cannot be computed because the solution is defined on a non-regularly spaced mesh.")
    else
        return sum(η .* v) * (x[2] - x[1])
    end
end

"""
$(TYPEDSIGNATURES)

Compute the excess of mass of a solved initial-value problem `pb` at a given time `T`.

Keyword argument `T` is optional, the last computed time is used by default.

"""
function energy(p::Problem; T = nothing)
    η, v, x = solution(p; T = T)
    if !(x[2:end] .- x[2] ≈ x[1:(end - 1)] .- x[1])
        @error("The energy cannot be computed because the solution is defined on a non-regularly spaced mesh.")
    end
    mesh = Mesh(x)
    U = p.model.mapto(Init(mesh, η, v))
    ∂ₓ = 1im * mesh.k
    p.model.f!(U);fftm = -U[1] ./ ∂ₓ;fftm[1] = 0
    m = real(ifft(fftm))
    return @. mesh.dx / 2 * ($sum(η^2) + $sum(v * m))
end

"""
$(TYPEDSIGNATURES)

Compute the difference of excess of mass of a solved initial-value problem `pb` between given time `T` and initial time.

Keyword argument `T` is optional, the last computed time is used by default.

If keyword argument `rel=true` (default is false), then compute the relative difference (with initial value as reference).

"""
function mass_diff(p::Problem; T = nothing, rel = false)
    η, v, x = solution(p; T = T)
    η0, v0, x0 = solution(p; T = 0)
    return if !(x[2:end] .- x[2] ≈ x[1:(end - 1)] .- x[1])
        @error("The excess of mass difference cannot be computed because the solution is defined on a non-regularly spaced mesh.")
    else
        if rel == false
            return sum(η - η0) * (x[2] - x[1])
        else
            return sum(η - η0) / sum(η0)
        end
    end
end

"""
$(TYPEDSIGNATURES)

Compute the difference of horizontal impulse of a solved initial-value problem `pb` between given time `T` and initial time.

Keyword argument `T` is optional, the last computed time is used by default.

If keyword argument `rel=true` (default is false), then compute the relative difference (with initial value as reference).

"""
function momentum_diff(p::Problem; T = nothing, rel = false)
    η, v, x = solution(p; T = T)
    η0, v0, x0 = solution(p; T = 0)
    return if !(x[2:end] .- x[2] ≈ x[1:(end - 1)] .- x[1])
        @error("The horizontal impulse difference cannot be computed because the solution is defined on a non-regularly spaced mesh.")
    else
        if rel == false
            return sum((η - η0) .* v + η0 .* (v - v0)) * (x[2] - x[1])
        else
            return sum((η - η0) .* v + η0 .* (v - v0)) / sum(η0 .* v0)
        end
    end
end

"""
$(TYPEDSIGNATURES)

Compute the difference of energy of a solved initial-value problem `pb` between given time `T` and initial time.

Keyword argument `T` is optional, the last computed time is used by default.

If keyword argument `rel=true` (default is false), then compute the relative difference (with initial value as reference).

"""
function energy_diff(p::Problem; T = nothing, rel = false)
    η, v, x = solution(p; T = T)
    η0, v0, x0 = solution(p; T = 0)
    if !(x[2:end] .- x[2] ≈ x[1:(end - 1)] .- x[1])
        @error("The energy difference cannot be computed because the solution is defined on a non-regularly spaced mesh.")
    end
    mesh = Mesh(x)
    U = p.model.mapto(Init(mesh, η, v))
    U0 = p.model.mapto(Init(mesh, η0, v0))
    ∂ₓ = 1im * mesh.k
    p.model.f!(U);fftm = -U[1] ./ ∂ₓ;fftm[1] = 0
    p.model.f!(U0);fftm0 = -U0[1] ./ ∂ₓ;fftm0[1] = 0
    m = real(ifft(fftm));    m0 = real(ifft(fftm0))
    if rel == false
        return @. mesh.dx / 2 * ($sum((η - η0) * η + η0 * (η - η0)) + $sum((v - v0) * m + v0 * (m - m0)))
    else
        return @. ($sum((η - η0) * η + η0 * (η - η0)) + $sum((v - v0) * m + v0 * (m - m0))) / ($sum(η0^2) + $sum(v0 * m0))
    end
end
