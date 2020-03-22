function SolitaryWaveWhitham(mesh :: Mesh, param :: NamedTuple, guess :: Vector{Float64}; iterative = false :: Bool, verbose = true :: Bool, max_iter = 50, tol = 1e-14, KdV = false)
        # A good guess for low velocities is
        #function sol(x,α)
        #	2*α*sech.(sqrt(3*2*α)/2*x).^2
        #end
        c = param.c
        ϵ = param.ϵ
        μ = param.μ

        k = mesh.k
        F₁ = sqrt.(tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k)))
        F₁[1]=1
        if KdV == true
                F₁ = 1 .-μ/6*k.^2
        end


        Π⅔ = abs.(k) .< maximum(k) * 2/3
        Dx       =  1im * mesh.k

        function proj( u :: Vector{Float64} )
            real.(ifft(Π⅔.*fft(u)))
        end

        function proj( u :: Vector{Float64}, v :: Vector{Float64} )
            u-0*(v'*u)*v/(norm(v,2)^2)
        end

        function F( u :: Vector{Float64} )
            -c*u+real.(ifft(F₁.*fft(u)))+3*ϵ/4*u.^2
        end

        if iterative == false
                k = mesh.k
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
                        du .= -proj( JacF(u) \ fu , dxu )
                else
                        du .= -proj( gmres( JacFfast(u) , fu ) , dxu )
                end
    		u .+= real.(du)
        end
        return (u,flag)
end
