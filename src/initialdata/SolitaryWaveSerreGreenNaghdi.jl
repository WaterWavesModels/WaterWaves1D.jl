export SolitaryWaveSerreGreenNaghdi
using LinearMaps,IterativeSolvers,LinearAlgebra

"""
    SolitaryWaveSerreGreenNaghdi(param; kwargs...)

Compute the Serre-Green-Naghdi solitary wave with prescribed velocity.

# Arguments
- `param :: NamedTuple`: parameters of the problem containing velocity `c` and dimensionless parameters `ϵ` and `μ`, and mesh size `L` and number of collocation points `N`;
- `x₀ :: Real`: (keyword, optional, default = 0) center of solitary wave.

# Return values
`(η,u,v)` with
- `η :: Vector{Float64}`: surface deformation;
- `u :: Vector{Float64}`: layer-averaged velocity;
- `v :: Vector{Float64}`: tangential velocity;
- `mesh :: Mesh`: mesh collocation points.

"""
function SolitaryWaveSerreGreenNaghdi(
                param :: NamedTuple;
                x₀ = 0 :: Real
                        )


        c = param.c
        ϵ = param.ϵ
        μ = param.μ

        mesh = Mesh(param)

        η = (c^2-1)/ϵ*sech.(sqrt(3*(c^2-1)/(c^2)/μ)/2*(mesh.x.-x₀)).^2

        h = 1 .+ ϵ*η
        u = c*η./h

	k = mesh.k
        Dx       =  1im * k

        DxF(v) = real.(ifft(Dx .* fft(v)))
	v = u - μ/3 ./h .* (DxF(h.^3 .*DxF(u)))

        return (η,u,v,mesh)

end
