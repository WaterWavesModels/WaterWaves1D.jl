function SolitaryWaveWhitham(mesh :: Mesh, param :: NamedTuple, guess :: Vector{Float64}; iterative = false :: Bool, verbose = true :: Bool, max_iter = 50, tol = 1e-14)
        # A good guess for low velocities is
        #function sol(x,α)
        #	2*α*sech.(sqrt(3*2*α)/2*x).^2
        #end
        c = param.c
        ϵ = param.ϵ
        μ = param.μ

        k = mesh.k
        F₁ 	= sqrt.(tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k)))
        F₁[1]=1

        function F( u :: Vector{Float64} )
            -c*u+real.(ifft(F₁.*fft(u)))+3*ϵ/4*u.^2
        end

        if iterative == false
                x = mesh.x
                x₀ = mesh.x[1]
                FFT = exp.(-1im*k*(x.-x₀)');
                IFFT = exp.(1im*k*(x.-x₀)')/length(x);
                M = real.(IFFT*(diagm( 0 => F₁)*FFT))
                function JacF( u₀ :: Vector{Float64} )
                    M+diagm(0 => 3*ϵ/2*u₀ .-c)
                end
        else
                function JacFfast( u₀ :: Vector{Float64} )
                        dF(u) = -c*u+real.(ifft(F₁.*fft(u)))+3*ϵ/2*u₀.*u
                        return LinearMap(dF, length(u₀); issymmetric=true, ismutating=false)
                end
        end


        iter = 0
        err = 1
        u = copy(guess)
        du = similar(F(u))
        fu = similar(F(u))
        for i in range(1, length=max_iter)
                fu .= F(u)
    	        err = norm(fu,2)/norm(u,2)
    		if err < tol
    			@info string("Converged : ",err,"\n")
    			break
    		elseif verbose == true
                        print(string("error at step ",i,": ",err,"\n"))
    		end
                if i == max_iter
                        @warn  "The algorithm did not converge"
                end
                if iterative == false
                        du .= -JacF(u) \ fu
                else
                        du .= -gmres(JacFfast(u) , fu)
                end
    		u .+= real.(du)
        end
        return u
end
