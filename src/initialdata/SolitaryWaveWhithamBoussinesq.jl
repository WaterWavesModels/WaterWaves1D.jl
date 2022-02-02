export SolitaryWaveWhithamBoussinesq, SolitaryWB

"""
    SolitaryWaveWhithamBoussinesq(param; kwargs)

Compute the Whitham-Boussinesq solitary wave with prescribed velocity.

# Argument
- `param :: NamedTuple`: parameters of the problem containing velocity `c` \
and dimensionless parameters `ϵ` and `μ`,
and mesh size `L` and number of collocation points `N`;

## Keywords (optional)
- `guess :: Vector{Real}`: initial guess for the surface deformation (if not provided, the exact formula for SGN is used);
- `x₀ :: Real`: center of solitary wave (if guess is not provided);
- `α :: Real`: determines the model used (typically `1` or `1/2`, default is 1);
- `Boussinesq`: if `true` (default is `false`), compute the standard Boussinesq system with parameters `a` (defaut `-1//3`), `b=d` (defaut `1//3`), and `c=0`);
- `iterative :: Bool`: inverts Jacobian through GMRES if `true`, LU decomposition if `false`;
- `verbose :: Bool`: prints numerical errors at each step if `true`;
- `max_iter :: Int`: maximum number of iterations of the Newton algorithm;
- `tol :: Real`: general tolerance (default is `1e-10`);
- `ktol :: Real`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `gtol :: Real`: relative tolerance of the GMRES algorithm (default is `1e-10`);
- `dealias :: Int`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `q :: Real`: Newton algorithm modified with
`u_{n+1}=q*u_{n+1}+(1-q)*u_n`
(default is `1`);
- `β :: Real`: adds `β` times spectral projection onto the Kernel to the Jacobian.β

# Return values
`(η,v,mesh)` with
- `η :: Vector{Float64}`: surface deformation;
- `v :: Vector{Float64}`: velocity (derivative of the trace of the velocity potential at the surface);
- `mesh :: Mesh`: mesh collocation points.

"""
function SolitaryWaveWhithamBoussinesq(
        param :: NamedTuple;
        guess = zeros(0) :: Vector{Float64},
        x₀ = 0 :: Real,
        α = 1 :: Real,
	a = -1//3, b = 1//3,
	Boussinesq = false :: Bool,
        iterative = false :: Bool,
        verbose = false :: Bool,
        max_iter = 20 :: Int,
        tol = 1e-10 :: Real,
        ktol = 0 :: Real,
        gtol = 1e-10 :: Real,
        dealias = 0 :: Int,
        q=1 :: Real,
        β=0 :: Real)

        c = param.c
        ϵ = param.ϵ
        μ = param.μ
        if abs(c)<1
                @error("The velocity must be greater than 1 (in absolute value).")
        end


        mesh=Mesh(param)
        if guess == zeros(0)
                # Using the exact formula for the SGN solitary wave as initial guess (if not provided)
                η = (c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*(mesh.x.-x₀)).^2
                h = 1 .+ param.ϵ*η
                u = c*η./h
                k = mesh.k
        	F₀ = sqrt(param.μ)*1im * k
        	DxF(v) = real.(ifft(F₀ .* fft(v)))
        	guess = u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))
        end

        k = mesh.k
	if Boussinesq==false
		F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
		F₁[1] 	= 1
		F₂ = F₁.^α
	else
		F₂ = 1 ./(1 .+μ*b*abs.(k).^2)
		F₁ 	= (1 .-μ*a*abs.(k).^2).*(F₂.^2)
	end
        Dx       =  1im * sqrt(μ)*k

        if dealias == 0
                Π = ones(size(k)) # no dealiasing (Π=Id)
        else
                K = (mesh.kmax-mesh.kmin)/(2+dealias)
                Π = abs.(k) .<= K # Dealiasing low-pass filter
        end
        krasny(k) = (abs.(k).>= ktol ).*k
        krasny!(k) = k[abs.(k).< ktol ].=0

        function filt( v :: Vector{Float64} )
                real.(ifft(krasny(Π.*fft(v))))
        end
        function filt( v :: Vector{Complex{Float64}} )
                ifft(krasny(Π.*fft(v)))
        end
        function Four1( v  )
                ifft(F₁.*fft(v))
        end
        function Four2( v  )
                ifft(F₂.*fft(v))
        end
        function FF( v )
                w = Four2(v)
                return real.(-c^2*v .+ Four1(v) .+ c*ϵ/2*w.^2 .+ ϵ*Four2(c*w.*v .- ϵ/2*w.^3))
        end
        function Fabs( v )
                w = Four2(v)
                return real.(c^2*abs.(v) .+ abs.(Four1(v)) .+ c*ϵ/2*w.^2 .+ ϵ*Four2(c*abs.(w.*v) .+ ϵ/2*abs.(w).^3))
        end

        if iterative == false
                k = mesh.k
                x = mesh.x
                x₀ = mesh.x[1]
                FFT = exp.(-1im*k*(x.-x₀)');
                IFFT = exp.(1im*k*(x.-x₀)')/length(x);
                Id = Diagonal( ones(size(x)));
                M₁ = IFFT*(Diagonal( F₁)*FFT)
                M₂ = IFFT*(Diagonal( F₂)*FFT)
                function JacF( v ,dxv )
                        M₀ = Diagonal( Four2(v))
                        Symmetric(real.(-c^2 .* Id .+ M₁ .+ c*ϵ.*(M₀*M₂ .+ M₂*M₀) .+ ϵ*M₂*(Diagonal(c*v -3*ϵ/2*Four2(v).^2))*M₂ .+ β*dxv*dxv'))
                end
        else
                function JacFfast( v , dxv )
                        dF(φ) = real.(-c^2*φ+Four1(φ)+c*ϵ*Four2(v).*Four2(φ)+ϵ*Four2(c*Four2(v).*φ+(c*v -3*ϵ/2*(Four2(v).^2)).*Four2(φ))+ β*dot(dxv,φ)*dxv)
                        return LinearMap(dF, length(v); issymmetric=true, ismutating=false)
                end
        end

        flag=0
        iter = 0
        u = filt(guess)
        du = similar(u)
        fu = similar(u)
        dxu = similar(u)
        for i in range(0, stop=max_iter)
                dxu .= real.(ifft(Dx.*fft(u)))
                dxu ./= norm(dxu,2)
                fu .= FF(u)
                relerr = norm(fu,Inf)/norm(Fabs(u),Inf)
                abserr = norm(fu,Inf)
                if relerr < tol
    			@info string("Converged : relative error ",relerr," in ",i," steps\n")
    			break
    		elseif verbose == true
                        print(string("absolute error at step ",i,": ",abserr,"\n"))
                        print(string("relative error at step ",i,": ",relerr,"\n"))

    		end
                if i == max_iter
                        @warn string("The algorithm did not converge after ",i," steps: final relative error is ",relerr,"\n")
                        break
                end

                if iterative == false
                        du .= JacF(u,dxu) \ fu
                else
                        du .= gmres( JacFfast(u,dxu) , fu ; verbose = verbose, reltol = gtol )
                end
    		u .-= q*filt(du)
        end
        return (c*u-ϵ/2*real.(Four2(u).^2),u,mesh)
end

"""
    SolitaryWB(param; kwargs)

Build the initial data associated with `SolitaryWaveWhithamBoussinesq(param; kwargs)`, of type `InitialData`,
to be used in initial-value problems `Problem(model, initial::InitialData, param)`.

---
    SolitaryWB(c; ϵ=1,μ=1,N=2^12,kwargs)
Build the initial data with velocity `c`, dimensionless parameters `ϵ` and `μ`, and number of collocation points `N`, \
and `kwargs` the other (optional) keyword arguments as above.
"""
struct SolitaryWB <: InitialData

	η
	v
	label :: String
	info  :: String

	function SolitaryWB(param::NamedTuple;
				guess = zeros(0) :: Vector{Float64},
			        x₀ = 0 :: Real,
			        α = 1 :: Real,
				a = -1//3, b = 1//3,
				Boussinesq = false :: Bool,
			        iterative = false :: Bool,
			        verbose = false :: Bool,
			        max_iter = 20 :: Int,
			        tol = 1e-10 :: Real,
			        ktol = 0 :: Real,
			        gtol = 1e-10 :: Real,
			        dealias = 0 :: Int,
			        q=1 :: Real,
			        β=0 :: Real)

                (η,v,mesh)=SolitaryWaveWhithamBoussinesq(param;
                                guess=guess,x₀=x₀,α=α,a=a,b=b,Boussinesq=Boussinesq,
                                iterative=iterative,verbose=verbose,max_iter=max_iter,
                                tol=tol,ktol=ktol,gtol=gtol,dealias=dealias,q=q,β=β)
		init = Init(mesh,η,v)
		if Boussinesq == false
        		model = "Whitham-Boussinesq"
		else
			model = "Boussinesq"
		end
                label = "$model solitary wave"
		info = "Solitary travelling wave for the $model model.\n\
		├─velocity c = $(param.c)\n\
		└─maximum h₀ = $(maximum(η)) (from rest state)."

		new( init.η,init.v,label,info  )
	end

	function SolitaryWB(
                        c::Real; ϵ=1::Real,μ=1::Real,N=2^12::Int,
			guess = zeros(0) :: Vector{Float64},
		        x₀ = 0 :: Real,
		        α = 1 :: Real,
			a = -1//3, b = 1//3,
			Boussinesq = false :: Bool,
		        iterative = false :: Bool,
		        verbose = false :: Bool,
		        max_iter = 20 :: Int,
		        tol = 1e-10 :: Real,
		        ktol = 0 :: Real,
		        gtol = 1e-10 :: Real,
		        dealias = 0 :: Int,
		        q=1 :: Real,
		        β=0 :: Real)

		L=200/sqrt(3*(c^2-1)/(c^2)/μ)/2
		xmin,xmax = x₀-L,x₀+L;
		param=(ϵ=ϵ,μ=μ,c=c,L=L,xmin=xmin,xmax=xmax,N=N)
		sol=SolitaryWB(param;
                        guess=guess,x₀=x₀,α=α,a=a,b=b,Boussinesq=Boussinesq,
                        iterative=iterative,verbose=verbose,max_iter=max_iter,
                        tol=tol,ktol=ktol,gtol=gtol,dealias=dealias,q=q,β=β)
		new( sol.η,sol.v,sol.label,sol.info  )
	end

end
