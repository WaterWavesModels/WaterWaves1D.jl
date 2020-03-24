function SolitaryWaveWhithamGreenNaghdi2c(mesh :: Mesh, param :: NamedTuple, guess :: Vector{Float64}; iterative = false :: Bool, verbose = true :: Bool, max_iter = 50, tol = 1e-14, q=1, GN = false)
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
        krasny(k) = (abs.(k).> 1e-12 ).*k
        krasny!(k) = k[abs.(k).< 1e-12 ].=0

        Dx       =  1im * mesh.k

        function dealias( v  )
                real.(ifft(krasny(Π⅔.*fft(v))))
        end
        #
        # function F( v , h , Fv )
        #         return -1/3 ./ (h.^2).* ifft(F₀.*fft( (h.^3).* Fv )) .+
        #                 ϵ/2 .*  (h.^2).*Fv.^2 .+
        #                 v .- 1/(c^2) .*v.*h .- ϵ/2 .* v.^2
        # end
        function F( v , h , Fv )
                return 1/3 *(v .-1)./ h .* ifft(F₀.*fft( (h.^3).* Fv )) .+
                        ϵ/2 .*  (h.^2).*Fv.^2 .+
                        v .- 1/(c^2) .*(h .-1) .- ϵ/2 .* v.^2
        end

        function Fabs( v, h, Fv )
                return abs.(1/3 ./ (h.^2).* ifft(F₀.*fft( (h.^3).* Fv))) .+
                        abs.(ϵ/2 .* ( h.*Fv ).^2) .+
                        abs.(v) .+ abs.(1/(c^2) .*v.*h) .+ abs.(ϵ/2 .* v.^2)
        end
        if iterative == false
                k = mesh.k
                x = mesh.x
                x₀ = mesh.x[1]
                FFT = exp.(-1im*k*(x.-x₀)');
                IFFT = exp.(1im*k*(x.-x₀)')/length(x);
                Id = Diagonal(ones(size(x)));
                M₀ = IFFT*(Diagonal( F₀)*FFT)
                M(v) = Diagonal( v )
                function JacF( v , h, h2, Fv )
                        -1/3 *M( 1 ./ h.^2 )* M₀ * M(h.^3)* ( M₀ .+
                                3*ϵ *  M( h .*Fv ) ) .+
                                 2*ϵ/3 * M( 1 ./ h .* ( ifft(F₀.*fft(h.^3 .*Fv ) ) ) ).+
                                 ϵ * M( h.^2 .* Fv ) * M₀ .+
                                 ϵ^2 * M( h.^3 .* Fv.^2 ) .+
                                 Id .- M(h.^2)./(c^2) .- ϵ* M( v )
                end
        else
                function JacFfast( v :: Vector{Float64} )
                        dF(φ) = -1/3 ./ (h(v).^2).* Four( (h(v).^3).* Four( φ ) ) .-
                                ϵ ./ (h(v).^2).* Four( (h(v).^4).* φ .* Four( v ) ) .+
                                 2*ϵ/3 .* φ ./ h(v).* Four( (h(v).^3).* Four( v ) ) .+
                                 ϵ .* (h(v).*Four(v)).*( h(v) .* Four( φ ) .+ ϵ.*φ.* (h(v).^2) .* Four(v) ) .+
                                 φ .- φ.* (h(v).^2)./(c^2) .- ϵ* v.* φ
                        return LinearMap(dF, length(v); issymmetric=true, ismutating=false)
                end
        end

        flag=0
        iter = 0
        err = 1
        u = fft(guess./(1 .+ ϵ*guess)) # initial guess for iteration
        du = similar(u)
        h = similar(u)
        h2 = similar(u)
        v = similar(u)
        fu = similar(u)
        dxu = similar(u)
        for i in range(1, length=max_iter)
                krasny!(u)
                v.=ifft(u)
                h.=1 ./(1 .-v)
                h2.=1 ./(1 .-v).^2
                dxu .= real.(ifft(Dx.*u))
                dxu ./= norm(dxu,2)
                fu .= F(v,h,ifft(F₀.*u))
                err = norm(fu,Inf)/norm(Fabs(v,h,ifft(F₀.*u)),Inf)
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
                        du .= JacF(v,h,h2,ifft(F₀.*u)) \ fu
                else
                        du .= gmres( JacFfast(u) , fu )
                end
    		u .-= fft(q*du)
                krasny!(u)
        end

        return (real.(ifft(u)./(1 .- ϵ*ifft(u))),c*real.(ifft(u)),flag)
end
