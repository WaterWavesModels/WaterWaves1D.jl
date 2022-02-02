export SolitaryWaveWhitham, SolitaryWhitham

"""
    SolitaryWaveWhitham(param; kwargs)

Compute the Whitham solitary wave with prescribed velocity.

# Argument
- `param :: NamedTuple`: parameters of the problem containing velocity `c` \
and dimensionless parameters `ϵ` and `μ`, \
and mesh size `L` and number of collocation points `N`;

## Keywords (optional)
- `guess :: Vector{Real}`: initial guess for the surface deformation (if not provided, the exact formula for KdV is used);
- `x₀ :: Real`: center of solitary wave (if guess is not provided);
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
- `α :: Real`: adds `α` times spectral projection onto the Kernel to the Jacobian;
- `KdV :: Bool`: if `true` computes the KdV (instead of Whitham) solitary wave.
# Return values
`(u,mesh)` with
- `u :: Vector{Float64}`: the solution;
- `mesh :: Mesh`: mesh collocation points.

"""
function SolitaryWaveWhitham(
                param :: NamedTuple;
                guess = zeros(0) :: Vector{Float64},
                x₀ = 0 :: Real,
                iterative = false :: Bool,
                verbose = false :: Bool,
                max_iter = 20 :: Int,
                tol = 1e-10 :: Real,
                gtol = 1e-10 :: Real,
                ktol = 0 :: Real,
                dealias = 0 :: Int,
                q=1 :: Real,
                α=0 :: Real,
                KdV = false :: Bool)


        mesh=Mesh(param)
        c = param.c
        ϵ = param.ϵ
        μ = param.μ
        if abs(c)<1
                @error("The velocity must be greater than 1 (in absolute value).")
        end


        if guess == zeros(0)
                # Using the exact formula for the KdV solitary wave as initial guess
                guess = 2*(c-1)/ϵ*sech.(sqrt(3*2*(c-1)/μ)/2*mesh.x.-x₀).^2
        end

        k  = mesh.k
        Dx =  1im * k
        if KdV == true
                F₁ = 1 .-μ/6*k.^2
        else
                F₁ = sqrt.(tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k)))
                F₁[1]=1
        end

        if dealias == 0
                Π = ones(size(k)) # no dealiasing (Π=Id)
        else
                K = (mesh.kmax-mesh.kmin)/(2+dealias)
                Π = abs.(k) .<= K # Dealiasing low-pass filter
        end
        krasny(k) = (abs.(k).> ktol ).*k
        krasny!(k) = k[abs.(k).< ktol ].=0

        function filter( v :: Vector{Float64} )
                real.(ifft(krasny(Π.*fft(v))))
        end
        function filter( v :: Vector{Complex{Float64}} )
                ifft(krasny(Π.*fft(v)))
        end


        function F( u :: Vector{Float64} )
            -c*u+real.(ifft(F₁.*fft(u)))+3*ϵ/4*u.^2
        end
        function Fabs( u :: Vector{Float64} )
            c*abs.(u)+abs.(ifft(F₁.*fft(u)))+3*ϵ/4*u.^2
        end

        if iterative == false
                k = mesh.k
                x = mesh.x
                x₀ = mesh.x[1]
                FFT = exp.(-1im*k*(x.-x₀)');
                IFFT = exp.(1im*k*(x.-x₀)')/length(x);
                M = real.(IFFT*(Diagonal( F₁)*FFT))
                function JacF( v , dxv  )
                    Symmetric(M+Diagonal(3*ϵ/2*v .-c)+ α*dxv*dxv')
                end
        else
                function JacFfast( v ,dxv )
                        dF(φ) = -c*φ+real.(ifft(F₁.*fft(φ)))+3*ϵ/2*v.*φ+ α*dot(dxv,φ)*dxv
                        return LinearMap(dF, length(v); issymmetric=true, ismutating=false)
                end
        end

        flag=0
        iter = 0
        u = filter(guess)
        du = similar(u)
        fu = similar(u)
        dxu = similar(u)
        for i in range(0, stop=max_iter)
                dxu .= real.(ifft(Dx.*fft(u)))
                dxu ./= norm(dxu,2)
                fu .= F(u)
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
    		u .-= q*filter(du)
        end
        return (u,mesh)
end


"""
    SolitaryWhitham(param; kwargs)

Build the initial data associated with `SolitaryWaveWhitham(param; kwargs)`, of type `InitialData`,
to be used in initial-value problems `Problem(model, initial::InitialData, param)`.

---
    SolitaryWhitham(c; ϵ=1,μ=1,N=2^12,kwargs)
Build the initial data with velocity `c`, dimensionless parameters `ϵ` and `μ`, and number of collocation points `N`, \
and `kwargs` the other (optional) keyword arguments as above.
"""
struct SolitaryWhitham <: InitialData

	η
	v
	label :: String
	info  :: String

	function SolitaryWhitham(param::NamedTuple;
			guess = zeros(0) :: Vector{Float64},
			x₀ = 0 :: Real,
			iterative = false :: Bool,
			verbose = false :: Bool,
			max_iter = 20 :: Int,
			tol = 1e-10 :: Real,
			gtol = 1e-10 :: Real,
			ktol = 0 :: Real,
			dealias = 0 :: Int,
			q=1 :: Real,
			α=0 :: Real,
			KdV = false :: Bool)
                (u,mesh)=SolitaryWaveWhitham(param;
                                guess=guess,x₀=x₀,KdV=KdV,
                                iterative=iterative,verbose=verbose,max_iter=max_iter,
                                tol=tol,ktol=ktol,gtol=gtol,dealias=dealias,q=q,α=α)
		ϵ=param.ϵ
		η=((ϵ*u/2 .+1).^2 .-1)/ϵ
		init = Init(mesh,η,u)
		if KdV == false
        		model = "Whitham"
		else
			model = "KdV"
		end
                label = "$model solitary wave"
		info = "Solitary travelling wave for the $model model.\n\
		├─velocity c = $(param.c)\n\
		└─maximum h₀ = $(maximum(η)) (from rest state)."

		new( init.η,init.v,label,info  )
	end

	function SolitaryWhitham(
                        c::Real; ϵ=1::Real,μ=1::Real,N=2^12::Int,
			guess = zeros(0) :: Vector{Float64},
	                x₀ = 0 :: Real,
	                iterative = false :: Bool,
	                verbose = false :: Bool,
	                max_iter = 20 :: Int,
	                tol = 1e-10 :: Real,
	                gtol = 1e-10 :: Real,
	                ktol = 0 :: Real,
	                dealias = 0 :: Int,
	                q=1 :: Real,
	                α=0 :: Real,
	                KdV = false :: Bool)
		L=200/sqrt(3*(c^2-1)/(c^2)/μ)/2
		xmin,xmax = x₀-L,x₀+L;
		param=(ϵ=ϵ,μ=μ,c=c,L=L,xmin=xmin,xmax=xmax,N=N)
		sol=SolitaryWhitham(param;
                                guess=guess,x₀=x₀,KdV=KdV,
                                iterative=iterative,verbose=verbose,max_iter=max_iter,
                                tol=tol,ktol=ktol,gtol=gtol,dealias=dealias,q=q,α=α)
		new( sol.η,sol.v,sol.label,sol.info  )
	end

end
