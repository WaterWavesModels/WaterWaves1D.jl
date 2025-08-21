export HBBMb

"""
    HBBMb(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the hyperbolized BBM system 

# Argument
`param` is of type `NamedTuple` and must contain
- relaxation parameter `ϵ` (nonlinearity) and `δ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `prep`: `∈{0,1,2}` and represent the level of preparation of the initial data (default is `1`);
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `stable`: consider the stable system if `stable=1`  (default is `1/2`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references;

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `HBBMb.f!` to be called in explicit time-integration solvers;
2. a function `HBBMb.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `HBBMb.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
	- `η` is the values of surface deformation at collocation points `x`;
	- `v` is the derivative of the trace of the velocity potential at `x`;
4. additionally, a handy function `HBBMb.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u,p,w)`, where
	- `u` corresponds to the layer-averaged horizontal velocity.
	- `p` corresponds to the relaxed (artificial) layer-averaged non-hydrostatic pressure;
	- `w` corresponds to the relaxed (artificial) layer-averaged vertical velocity.

"""
mutable struct HBBMb <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info    :: String

    function HBBMb(param::NamedTuple; 
								prep = 1,
								mesh = Mesh(param),
								dealias	= 0,
								ktol	= 0,
								stable	= 1/2,
								label	= nothing
								)
		# Set up
		δ 	= param.δ
		ϵ 	= param.ϵ

		if isnothing(label)
				label = "hyperbolized BBM system"
		end

		# Print information
		info = "$label model.\n"
		info *= "├─Relaxation parameter ϵ=$ϵ.\n"
		info *= "├─Shallowness parameter δ=$δ.\n"
		info *= "├─Initial data prepared of order $prep.\n"
		if dealias == 0
			info *= "├─No dealiasing. "
		else
			info *= "├─Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k = mesh.k
		x 	= mesh.x
		∂ₓ	=  1im * k

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		u = zeros(Complex{Float64}, mesh.N)
		fftu, fftv, fftw, B = (similar(u),).*ones(4)

		# Evolution equations are ∂t U = f(U)
		function f!(U;ϵ=ϵ)
				fftu .= U[1]; fftv .= U[2]; fftw .= U[3] ; B .= U[4] ;
				u .= ifft(fftu)
	
				U[1] .= -∂ₓ.* (fftu .+ 1/2*Π⅔.*fft( u.^2 ) ) -  δ^2*fft(B.^(4*stable) .* ifft(∂ₓ.* fftv))
				U[2] .=  (fftw .- ∂ₓ.* fftu)/ϵ^2
				U[3] .= -fft(B.^(4*(1-stable)) .*ifft(fftv))
				U[4] .= zero(B)
				for u in U u[ abs.(u).< ktol ].=0 end
		end

		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			fftu .= fft(data.v(x))
			
			if prep == 0 
				fftw .=  fft(zero(x))
			else
				fftw .= ∂ₓ.*fftu
			end

			if prep <= 1
				fftv .= fft(zero(x))
			else
				u .= ifft(fftu)
				fftv .= (∂ₓ.^2).* (fftu .+ 1/2*Π⅔.*fft( u.^2 ))./ (1 .+ (δ*k).^2)				
			end
			
			U = [Π⅔.*fftu, Π⅔.*fftv, Π⅔.*fftw, 1 .- data.η(x)]

			for u in U u[ abs.(u).< ktol ].=0 end
			return U
		end

		# Reconstruct variables from raw data
		# Return `(η,v,x)`, where
		# - `u` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			fftu .= U[1];
			fftv .= U[2];
			fftw .= U[3]; 
			B.= U[4];
			real(ifft(fftu)),real(B),mesh.x
		end

        function mapfrofull(U)
			fftu .= U[1]
			fftv .= U[2]
			fftw .= U[3]; 
			B.= U[4];

			real(ifft(fftu)),real(ifft(fftv)),real(ifft(fftw)),real(B),mesh.x
		end

		new(label, f!, mapto, mapfro, mapfrofull, info )
    end


end
