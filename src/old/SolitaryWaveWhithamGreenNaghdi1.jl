function SolitaryWaveWhithamGreenNaghdi1(mesh :: Mesh, param :: NamedTuple, guess :: Vector{Float64}; iterative = false :: Bool, verbose = true :: Bool, max_iter = 50, tol = 1e-14, q=1, GN = false)
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
        F₀       = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(k)

        if GN == true
                F₀ = 1im * k ./ (1 .+ μ/3 * k.^2 ).^(1/4)
        end

        Π⅔ = abs.(k) .< maximum(k) * 4/3
        krasny(k) = (abs.(k).> 1e-32 ).*k
        krasny!(k) = k[abs.(k).< 1e-32 ].=0

        Dx       =  1im * mesh.k

        function dealias( v :: Vector{Float64} )
                real.(ifft(krasny(Π⅔.*fft(v))))
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
        function Four( v )
                real.(ifft(F₀.*fft(v)))
        end

        function E( η :: Vector{Float64} , u :: Vector{Float64} )
                return η.^2 .+ h(η).*(u.^2) .+ μ/3 .* (h(u).^3).* (Four( u ).^2)
        end

        function F( v :: Vector{Float64} )
                return -1/3 ./ (h(v).^2).* Four( (h(v).^3).* Four(η(v))) .+
                        ϵ/2 .* (h(v).*Four(η(v))).^2 .+
                        η(v) .- 1/(c^2) .*v .- ϵ/2 .* η(v).^2
        end

        function Fabs( v :: Vector{Float64} )
                return abs.(1/3 ./ (h(v).^2).* Four( (h(v).^3).* Four(η(v)))) .+
                        abs.(ϵ/2 .* (h(v).*Four(η(v))).^2 ) .+
                        abs.(η(v)) .- 1/(c^2) .*abs.(v) .+ ϵ/2 .* abs.(η(v).^2)
        end

        if iterative == false
                k = mesh.k
                x = mesh.x
                x₀ = mesh.x[1]
                FFT = exp.(-1im*k*(x.-x₀)');
                IFFT = exp.(1im*k*(x.-x₀)')/length(x);
                Id = Diagonal(ones(size(x)));
                M₀ = real.(IFFT*(Diagonal( F₀)*FFT))
                M(v) = Diagonal( v )
                function JacF( v :: Vector{Float64} )
                        -1/3 *M(1 ./ (h(v).^2))* M₀ * M(h(v).^3)* M₀ * M( 1 ./ (h(v).^2) ) .-
                                ϵ * M(1 ./ (h(v).^2)) * M₀ * M( (h(v).^2) .* (M₀*η(v)) ) .+
                                 2*ϵ/3 * M( 1 ./ (h(v).^3).* ( M₀*((h(v).^3).* (M₀*η(v)) ) ) ) .+
                                 ϵ * M( (h(v).^2) .* (M₀*η(v)) ) * M₀ * M( 1 ./ (h(v).^2) ) .+
                                 ϵ^2 * M( h(v).* (M₀ * η(v)).^2 ) .+
                                 M( 1 ./ (h(v).^2) ) .- Id./(c^2) .- ϵ* M( η(v) ./ (h(v).^2) )
                end
        else
                function JacFfast( v :: Vector{Float64} )
                        dF(φ) = -1/3 ./ (h(v).^2).* Four( (h(v).^3).* Four( φ ./ (h(v).^2) ) ) .-
                                ϵ ./ (h(v).^2).* Four( (h(v).^2).* φ .* Four( η(v) ) ) .+
                                 2*ϵ/3 .* φ ./ (h(v).^3).* Four( (h(v).^3).* Four( η(v) ) ) .+
                                 ϵ .* (h(v).*Four(η(v))).*( h(v) .* Four( φ ./ (h(v).^2) ) .+ ϵ.*φ .* Four(η(v)) ) .+
                                 φ ./ (h(v).^2) .- φ./(c^2) .- ϵ* η(v).* φ ./ (h(v).^2)
                        return LinearMap(dF, length(v); issymmetric=true, ismutating=false)
                end
        end

        flag=0
        iter = 0
        err = 1
        u = copy(dealias(guess))
        du = similar(F(u))
        fu = similar(F(u))
        dxu = similar(u)
        for i in range(1, length=max_iter)
                dxu .= real.(ifft(Dx.*fft(u)))
                dxu ./= norm(dxu,2)
                fu .= F(dealias(u))
    	        err = norm(fu,2)/sqrt(norm(E(u,c*η(u))/c^2,1))
                err = norm(fu,Inf)/norm(Fabs(u),Inf)
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
                        du .=  Symmetric(JacF(u).-dxu*dxu') \ fu
                else
                        du .=  gmres( JacFfast(u) , fu )
                end
    		u .-= dealias(q*du)
        end

        return (u,c*η(u),flag)
end
