function SolitaryWaveWhithamBoussinesq(mesh :: Mesh, param :: NamedTuple, guess :: Vector{Float64}; iterative = false :: Bool, verbose = true :: Bool, max_iter = 50, tol = 1e-14, KdV = false, q=1)
        # A good guess for low velocities is
        #function sol(x,α)
        #	2*α*sech.(sqrt(3*2*α)/2*x).^2
        #end
        c = param.c
        ϵ = param.ϵ
        μ = param.μ

        k = mesh.k
        F₁ 	= tanh.(sqrt(μ)*abs.(mesh.k))./(sqrt(μ)*abs.(mesh.k))
        F₁[1] 	= 1
        F₂ = F₁.^(param.α)

        Π⅔ = abs.(k) .< maximum(k) * 2/3
        Dx       =  1im * mesh.k

        function proj( u :: Vector{Float64} )
            real.(ifft(Π⅔.*fft(u)))
        end

        function proj( u :: Vector{Float64}, v :: Vector{Float64} )
            u-0*(v'*u)*v/(norm(v,2)^2)
        end

        function vtou( v :: Vector{Float64} )
                real.(ifft(F₂.*fft(v)))
        end
        function F( v :: Vector{Float64} )
                u = vtou(v)
                return -c^2*v+real.(ifft(F₁.*fft(v)))+c*ϵ/2*u.^2+ϵ*real.(ifft(F₂.*fft(c*u.*v-ϵ/2*u.^3)))
        end

        if iterative == false
                k = mesh.k
                x = mesh.x
                x₀ = mesh.x[1]
                FFT = exp.(-1im*k*(x.-x₀)');
                IFFT = exp.(1im*k*(x.-x₀)')/length(x);
                Id = diagm( 0 => ones(size(x)));
                M₁ = real.(IFFT*(diagm( 0 => F₁)*FFT))
                M₂ = real.(IFFT*(diagm( 0 => F₂)*FFT))
                function JacF( v₀ :: Vector{Float64} )
                        M₀ = diagm( 0 => vtou(v₀))
                        -c^2 * Id + M₁ + c*ϵ*(M₀*M₂ + M₂*M₀)+ ϵ*M₂*(c*diagm( 0 => v₀) -3*ϵ/2*M₀*M₀)*M₂
                end
        else
                function JacFfast( v₀ :: Vector{Float64} )
                        u₀ = vtou(v₀)
                        dF(v) = -c^2*v+real.(ifft(F₁.*fft(v)))+c*ϵ*u₀.*real.(ifft(F₂.*fft(v)))+ϵ*real.(ifft(F₂.*fft(c*u₀.*v+c*v₀.*real.(ifft(F₂.*fft(v))) -3*ϵ/2*(u₀.^2).*real.(ifft(F₂.*fft(v))))))
                        return LinearMap(dF, length(v₀); issymmetric=true, ismutating=false)
                end
        end

        flag=0
        iter = 0
        err = 1
        u = copy(guess)
        du = similar(F(u))
        fu = similar(F(u))
        dxu = similar(u)
        for i in range(1, length=max_iter)
                dxu .= real.(ifft(Dx.*fft(u)))
                fu .= proj(F(u),dxu)
    	        err = norm(fu,2)/norm(u,2)
    		if err < tol
    			@info string("Converged : ",err,"\n")
    			break
    		elseif verbose == true
                        print(string("error at step ",i,": ",err,"\n"))
    		end
                if i == max_iter
                        flag=1
                        @warn  "The algorithm did not converge"
                end
                if iterative == false
                        du .= proj( JacF(u) \ fu , dxu )
                else
                        du .= proj( gmres( JacFfast(u) , fu ) , dxu )
                end
    		u .-= q*du
        end
        η = c*u-ϵ/2*vtou(u).^2
        return (η,u,flag)
end
