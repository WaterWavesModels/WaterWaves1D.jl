function SolitaryWaveWhithamGreenNaghdi(mesh :: Mesh, param :: NamedTuple, guess :: Vector{Float64}; iterative = false :: Bool, verbose = true :: Bool, max_iter = 50, tol = 1e-14, q=1)
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
        F₀       = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(mesh.k)


        Π⅔ = abs.(k) .< maximum(k) * 2/3
        Dx       =  1im * mesh.k

        function proj( u :: Vector{Float64} )
            real.(ifft(Π⅔.*fft(u)))
        end

        function proj( u :: Vector{Float64}, v :: Vector{Float64} )
            u-0*(v'*u)*v/(norm(v,2)^2)
        end

        function h( v :: Vector{Float64} )
                1 .+ ϵ*v
        end
        function η( v :: Vector{Float64} )
                v./h(v)
        end
        function Four( v :: Vector{Float64} )
                real.(ifft(F₀.*fft(v)))
        end
        function F( v :: Vector{Float64} )
                return -1/3 ./ (h(v).^2).* Four( (h(v).^3).* Four(η(v))) .+ ϵ/2 .* (h(v).*Four(η(v))).^2 .+ η(v) .- 1/(c^2) .*v .- ϵ/2 .* η(v).^2
        end

        if iterative == false
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
                function JacFfast( v :: Vector{Float64} )
                        dF(φ) = -1/3 ./ (h(v).^2).* Four( (h(v).^3).* Four( φ ./ (h(v).^2) ) ) .-
                                ϵ ./ (h(v).^2).* Four( (h(v).^2).* φ .* Four( η(v) ) ) .+
                                 2*ϵ/3 .* φ ./ (h(v).^3).* Four( (h(v).^3).* Four( η(v) ) ) .+
                                 ϵ .* (h(v).*Four(η(v))).*( h(v) .* Four( φ ./ (h(v).^2) ) + ϵ.*φ .* Four(η(v)) ) .+
                                 φ ./ (h(v).^2) .- φ./(c^2) .- ϵ* η(v).* φ ./ (h(v).^2)
                        return LinearMap(dF, length(v); issymmetric=false, ismutating=false)
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

        return (u,c*η(u),flag)
end
