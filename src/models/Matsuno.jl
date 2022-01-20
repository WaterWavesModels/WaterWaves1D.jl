export Matsuno_fast,Matsuno

"""
	Matsuno_fast(param;dealias,verbose)

Same as `Matsuno`, but faster.
"""
mutable struct Matsuno_fast <: AbstractModel

	f!		:: Function
	mapto	:: Function
	mapfro	:: Function

    function Matsuno_fast(param::NamedTuple; dealias=true, verbose=true)
		if verbose
			@info "Build the Matsuno model."
			@warn "The velocity is consistent with the tangential velocity used for water waves
			only when they are null."
		end


		ϵ        = param.ϵ
        mesh     = Mesh(param)
        x        = mesh.x
        Γ        = abs.(mesh.k)
        Dx       =  1im * mesh.k        # Differentiation
        H        = -1im * sign.(mesh.k) # Hilbert transform
		if dealias == true || dealias == 1
			Π⅔    = Γ .< mesh.kmax * 2/3     # Dealiasing low-pass filter
		else
			Π⅔    = zero(Γ) .+ 1     # Dealiasing low-pass filter
		end

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₀ = zeros(Complex{Float64}, mesh.N)
        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)

        Px  = plan_fft(hnew; flags = FFTW.MEASURE)

		function f!(U)

		    for i in eachindex(hnew)
		        hnew[i] = Γ[i] * U[i,1]
		    end

		    ldiv!(unew, Px, hnew )

		    for i in eachindex(hnew)
		        hnew[i] = Dx[i] * U[i,1]
		    end

		    ldiv!(I₁, Px, hnew)

		    unew  .*= I₁

		    mul!(I₁, Px, unew)

		    I₁  .*= ϵ .* Π⅔
		    I₁  .-= hnew

		    ldiv!(hnew, Px, view(U,:,1))
		    ldiv!(unew, Px, view(U,:,2))

		    I₂    .= hnew .* unew

		    mul!(I₃, Px, I₂)

		    I₃    .*= Dx

		    for i in eachindex(H)
		        U[i,1]  = H[i] * U[i,2]
		        I₀[i] = Γ[i] * U[i,2]
		    end

		    ldiv!(I₂, Px, I₀)

		    I₂    .*= hnew

		    mul!(hnew, Px, I₂)

		    hnew  .*= H
		    I₃    .+= hnew
		    I₃    .*= ϵ .* Π⅔

		    for i in eachindex(I₃)
		        U[i,1] -= I₃[i]
		    end

		    I₃    .=  unew.^2

		    mul!(unew, Px, I₃)

		    unew  .*= Dx
		    unew  .*= ϵ/2 .* Π⅔
		    I₁    .-= unew

		    for i in eachindex(I₁)
		        U[i,2] =  I₁[i]
		    end

		end

		# This function is correct only when the initial data for v is zero.
		function mapto(data::InitialData)

		    [Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]

		end

		# This function is correct only when the initial data for v is zero.
		function mapfro(U)
		    real(ifft(view(U,:,1))),real(ifft(view(U,:,2)))
		end

		new(f!, mapto, mapfro )
    end
end

"""
	Matsuno(param;dealias,verbose)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the quadratic deep-water model proposed by [Matsuno](https://doi.org/10.1103/PhysRevLett.69.609).

# Arguments
`param` is of type `NamedTuple` and must contain
- the dimensionless parameters `ϵ` (nonlinearity);
- numerical parameters to construct the mesh of collocation points as `mesh = Mesh(param)`.

## Optional keyword arguments
- `dealias`: dealiasing with `1/3` Orlicz rule if `true` (by default) or no dealiasing if `false`;
- `verbose`: prints information if `true` (default is `true`).


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `DeepQuadratic.f!` to be called in explicit time-integration solvers;
2. a function `DeepQuadratic.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `DeepQuadratic.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v)`, where
    - `η` is the surface deformation;
    - `v` is a velocity variable which is *not* the derivative of the trace of the velocity potential (if not null).

"""
mutable struct Matsuno <: AbstractModel

	f!		:: Function
	mapto	:: Function
	mapfro	:: Function


    function Matsuno(param::NamedTuple ; dealias=true, verbose=true)
		if verbose
			@info "Build the Matsuno model."
			@warn "The velocity is consistent with the tangential velocity used for water waves
			only when they are null."
		end

		ϵ 	= param.ϵ
		mesh = Mesh(param)
		x   = mesh.x
        Γ 	= abs.(mesh.k)
    	∂ₓ	=  1im * mesh.k            # Differentiation
        H 	= -1im * sign.(mesh.k)     # Hilbert transform
		if dealias == true || dealias == 1
			Π⅔    = Γ .< mesh.kmax * 2/3     # Dealiasing low-pass filter
		else
			Π⅔    = zero(Γ) .+ 1     # Dealiasing low-pass filter
		end

        hnew = zeros(Complex{Float64}, mesh.N)
        unew = zeros(Complex{Float64}, mesh.N)

        I₁ = zeros(Complex{Float64}, mesh.N)
        I₂ = zeros(Complex{Float64}, mesh.N)
        I₃ = zeros(Complex{Float64}, mesh.N)


		function f!(U)

		   hnew .= ifft(U[:,1])
		   unew .= ifft(U[:,2])
		   I₃ .= fft(ifft((∂ₓ).*U[:,1]).*ifft((Γ).*U[:,1]))
		   I₁ .= H.*U[:,2].-ϵ*Π⅔.*(H.*fft(hnew.*ifft(Γ.*U[:,2])).+∂ₓ.*fft(hnew.*unew))
		   I₂ .= -(∂ₓ.*U[:,1])+ϵ*Π⅔.*(I₃-∂ₓ.*fft(unew.^2)/2)
		   #
		   U[:,1] .= I₁
		   U[:,2] .= I₂

		end

		# This function is correct only when the initial data for v is zero.
		function mapto(data::InitialData)

			[Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]

		end

		# This function is correct only when the initial data for v is zero.
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2]))
		end

		new(f!, mapto, mapfro )
    end
end
