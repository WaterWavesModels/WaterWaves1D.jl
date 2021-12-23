export Random
using FFTW
"""
    Random(;args)

Randomly generated initial data, based on provided (optional arguments) :
- `L` is the typical wavelength (default is `L=1`)
- `s` is the (real) Sobolev index regularity (default is `s=∞`)
- `λ` is the length of spatial localization (default is none)

The initial data `(η,v)` are generated through randomly chosen Fourier coefficients,
multiplied with weigth `w=min(1,|k/L|^(s+1/2))` or `w=min(1,exp(-|k/L|)` if `s=∞`.
If `λ` is provided, the function in spatial variables is multiplied by `exp(-|x/λ|^2)`,
and in any case normalized to have maximum absolute value 1.


"""
struct Random <: InitialData

    η
    v

    function Random(;L=1,s=Inf,λ=2)


        function generate( x )


            k = Mesh( x ).k
            if s == Inf
                w = 1 ./ exp.(log(10)*abs.(k*L/(2*π)))
            else
                w = 1 ./( 1 .+ 9*abs.(k*L/(2*π)).^(s+1/2) )
            end
            if λ == nothing
                φ = zero(x).+1
            else
                φ = exp.(-abs.((x/λ).^2))
            end

            θ = 2*π*rand(Float64,length(x))
            r = rand(Float64,length(x))
            Fourier = r.*exp.(-1im.*θ).*w
            Physic = real.(ifft(Fourier)).*φ
            return Physic./maximum(abs.(Physic))
        end
        η( x ) = generate( x )
        v( x ) = generate( x )


    	new( η,v )

    end

end
