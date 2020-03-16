export SolitaryWaveWhitham

"""
    SolitaryWaveWhitham(mesh,param)
    mesh is of type Mesh, constrcucted with mesh=Mesh(p) where p is a set of parameters containing L and N
    param should contain a value c (for the velocity), μ (shallow water parameter)
    guess is the initial guess

"""
struct SolitaryWaveWhitham <: InitialData

    η
    v

    function SolitaryWaveWhitham(mesh :: Mesh, p :: NamedTuple, guess :: Vector{Float64})
        ϵ = param.ϵ
        μ = mesh.μ
        k = mesh.k
        F₁ 	= tanh.(sqrt(μ)*abs.(mesh.k))./(sqrt(μ)*abs.(mesh.k))
        function F( u :: Vector{Float64} )
            -c*u+ifft(F₁.*fft(u))+3*ϵ/4*u.^2
        end
        function JacF( u₀ :: Vector{Float64}, u :: Vector{Float64} )
            -c*u+ifft(F₁.*fft(u))+3*ϵ/2*u₀.*u
        end

        max_iter = 1000
        tol = 1e-12
        iter = 0
        err = 1
        u = guess
        while iter>max_iter & err > tol
            du = JacF(u) \ F(u)
            u .+= du
            err = norm(F(u),2)/norm(u,2)
            iter .+=1
        end
        η = u
        v = zeros(length(mesh.x))
        Init(mesh, η, v)


        #new( η , v )

    end

end
