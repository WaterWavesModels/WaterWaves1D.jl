function SolitaryWaveWhithamGreenNaghdi1b(mesh :: Mesh, param :: NamedTuple, guess :: Vector{Float64}; iterative = false :: Bool, verbose = true :: Bool, max_iter = 50, tol = 1e-14, q=1, GN = false)
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


        function F( v , hv, Fv, F2v )
                return -1/3 ./ hv.^2 .* F2v .+
                        ϵ/2 .* (hv.*Fv).^2 .+
                        v./hv .- 1/(c^2) .*v .- ϵ/2 .* (v./hv).^2
        end

        function Fabs( v , hv, Fv, F2v )
                return abs.(1/3 ./ hv.^2 .* F2v) .+
                        abs.(ϵ/2 .* (hv.*Fv).^2 ) .+
                        abs.(v./hv) .- 1/(c^2) .*abs.(v) .+ ϵ/2 .* abs.((v./hv).^2)
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
                function JacF( v , hv, Fv , F2v )
                        -1/3 *M(1 ./ hv.^2)* M₀ * M(hv.^3)* ( M₀ * M( 1 ./ hv.^2 ) .+
                                3* ϵ * M( Fv ./hv ) ) .+
                                 ϵ * M( hv.^2 .* Fv ) * M₀ * M( 1 ./ hv.^2 ) .+
                                 M( 1 ./ hv.^3 .+ 2*ϵ/3 ./ hv.^3 .* F2v .+ ϵ^2 * hv.* Fv.^2) .- Id./(c^2)
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
        du = similar(u)
        fu = similar(u)
        Fu = similar(u)
        F2u = similar(u)
        hu = similar(u)
        dxu = similar(u)
        for i in range(1, length=max_iter)
                u.=dealias(u)
                dxu .= real.(ifft(Dx.*fft(u)))
                dxu ./= norm(dxu,2)
                hu .= 1 .+ϵ*u
                Fu .= real.(ifft(F₀.*fft(u./hu)))
                F2u .= real.(ifft(F₀.*fft( hu.^3 .* Fu ) ) )
                fu .= F(u,hu,Fu,F2u)
                err = norm(fu,Inf)/norm(Fabs(u,hu,Fu,F2u),Inf)
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
                        du .=  Symmetric(JacF(u, hu, Fu,F2u).-dxu*dxu') \ fu
                else
                        du .=  gmres( JacFfast(u) , fu )
                end
    		u .-= dealias(q*du)
        end

        return (u,c*u./(1 .+ϵ*u),flag)
end
