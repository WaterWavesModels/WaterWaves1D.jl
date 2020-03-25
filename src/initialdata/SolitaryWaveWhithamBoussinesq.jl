export SolitaryWaveWhithamBoussinesq
"""
    `SolitaryWaveWhithamBoussinesq(mesh, param, guess; kwargs...)`

Computes the Whitham-Boussinesq solitary wave with prescribed velocity.

# Arguments
- `mesh :: Mesh`: parameters of the numerical grid, e.g constructed through Mesh(L,N);
- `param :: NamedTuple`: parameters of the problem containing velocity `c` and dimensionless parameters `ϵ` and `μ`;
- `guess :: Vector{Float64}`: initial guess for the surface deformation.
## Keywords
- `model :: Real`: determines the model used (typically `1` or `1/2`, default is 1);
- `iterative :: Bool`: inverts Jacobian through GMRES if `true`, LU decomposition if `false`;
- `verbose :: Bool`: prints numerical errors at each step if `true`;
- `max_iter :: Int`: maximum number of iterations of the Newton algorithm;
- `tol :: Real`: general tolerance (default is `1e-10`);
- `ktol :: Real`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `gtol :: Real`: relative tolerance of the GMRES algorithm;
- `dealias :: Int`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `q :: Real`: Newton algorithm modified with
`u_{n+1}=q*u_{n+1}+(1-q)*u_n`
(default is `1`);
- `α :: Real`: adds `α` times spectral projection onto the Kernel to the Jacobian.
# Return values
`(η,u) :: Tuple{Vector{Float64},Vector{Float64}}`
with
- `η`: surface deformation;
- `u`: velocity.
"""
function SolitaryWaveWhithamBoussinesq(mesh :: Mesh,
        param :: NamedTuple,
        guess :: Vector{Float64};
        model = 1 :: Real,
        iterative = false :: Bool,
        verbose = true :: Bool,
        max_iter = 20 :: Int,
        tol = 1e-10 :: Real,
        ktol = 0 :: Real,
        gtol = 1e-10 :: Real,
        dealias = 0 :: Int,
        q=1 :: Real,
        α=0 :: Real)
        # A good guess for low velocities is
        #function sol(x,α)
        #	2*α*sech.(sqrt(3*2*α)/2*x).^2
        #end
        c = param.c
        ϵ = param.ϵ
        μ = param.μ

        k = mesh.k
        F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
        F₁[1] 	= 1
        F₂ = F₁.^(model)
        Dx       =  1im * sqrt(μ)*k

        if dealias == 0
                Π = ones(size(k))
        else
                Π = abs.(k) .< maximum(k) * (1-dealias/(dealias+2))
        end
        krasny(k) = (abs.(k).>= ktol ).*k
        krasny!(k) = k[abs.(k).< ktol ].=0

        function filter( v :: Vector{Float64} )
                real.(ifft(krasny(Π.*fft(v))))
        end
        function filter( v :: Vector{Complex{Float64}} )
                ifft(krasny(Π.*fft(v)))
        end
        function Four1( v  )
                ifft(F₁.*fft(v))
        end
        function Four2( v  )
                ifft(F₂.*fft(v))
        end
        function F( v )
                w = Four2(v)
                return real.(-c^2*v .+ Four1(v) .+ c*ϵ/2*w.^2 .+ ϵ*Four2(c*w.*v .- ϵ/2*w.^3))
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
                        Symmetric(real.(-c^2 .* Id .+ M₁ .+ c*ϵ.*(M₀*M₂ .+ M₂*M₀) .+ ϵ*M₂*(Diagonal(c*v -3*ϵ/2*Four2(v).^2))*M₂ .+ α*dxv*dxv'))
                end
        else
                function JacFfast( v , dxv )
                        dF(φ) = real.(-c^2*φ+Four1(φ)+c*ϵ*Four2(v).*Four2(φ)+ϵ*Four2(c*Four2(v).*φ+(c*v -3*ϵ/2*(Four2(v).^2)).*Four2(φ))+ α*dot(dxv,φ)*dxv)
                        return LinearMap(dF, length(v); issymmetric=true, ismutating=false)
                end
        end

        flag=0
        iter = 0
        err = 1
        u = filter(guess)
        du = similar(u)
        fu = similar(u)
        dxu = similar(u)
        for i in range(1, length=max_iter)
                dxu .= real.(ifft(Dx.*fft(u)))
                dxu ./= norm(dxu,2)
                fu .= F(u)
    	        err = norm(fu,Inf)
    		if err < tol
    			@info string("Converged : ",err,"\n")
    			break
    		elseif verbose == true
                        print(string("error at step ",i,": ",err,"\n"))
    		end
                if i == max_iter
                        flag=1
                        @warn  string("The algorithm did not converge : ",err,"\n")
                end
                if iterative == false
                        du .= JacF(u,dxu) \ fu
                else
                        du .= gmres( JacFfast(u,dxu) , fu ; verbose = verbose, tol = gtol )
                end
    		u .-= q*filter(du)
        end
        return (c*u-ϵ/2*real.(Four2(u).^2),u)
end
