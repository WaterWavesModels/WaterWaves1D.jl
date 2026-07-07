export SolitaryWaveSerreGreenNaghdi, SolitarySGN

"""
    SolitaryWaveSerreGreenNaghdi(param; x₀=0)

Compute the Serre-Green-Naghdi solitary wave with prescribed velocity.

# Arguments
- `param :: NamedTuple`: parameters of the problem containing velocity `c` \
and dimensionless parameters `ϵ` and `μ`, \
and mesh size `L` and number of collocation points `N`;
- `x₀ :: Real`: (keyword, optional, default = 0) center of solitary wave.

# Return values
`(η,u,v,mesh)` with
- `η :: Vector{Float64}`: surface deformation;
- `u :: Vector{Float64}`: layer-averaged velocity;
- `v :: Vector{Float64}`: derivative of the trace of the velocity potential at the surface;
- `mesh :: Mesh`: mesh collocation points.

"""
function SolitaryWaveSerreGreenNaghdi(
        param::NamedTuple;
        x₀ = 0::Real
    )


    c = param.c
    ϵ = param.ϵ
    μ = param.μ
    if abs(c) < 1
        @error("The velocity must be greater than 1 (in absolute value).")
    end


    mesh = Mesh(param)

    η = (c^2 - 1) / ϵ * sech.(sqrt(3 * (c^2 - 1) / (c^2) / μ) / 2 * (mesh.x .- x₀)) .^ 2

    h = 1 .+ ϵ * η
    u = c * η ./ h

    k = mesh.k
    Dx = 1im * k

    DxF(v) = real.(ifft(Dx .* fft(v)))
    v = u - μ / 3 ./ h .* (DxF(h .^ 3 .* DxF(u)))

    return (η, u, v, mesh)

end

"""
    SolitarySGN(param; x₀=0)

Build the initial data associated with `SolitaryWaveSerreGreenNaghdi(param; x₀=0)`, of type `InitialData`,
to be used in initial-value problems `Problem(model, initial::InitialData, param)`.

---
	SolitarySGN(c; ϵ=1,μ=1,x₀=0,N=2^12)

Build the initial data with velocity `c`, center `x₀`, dimensionless parameters `ϵ` and `μ`, and number of collocation points `N`.
"""
struct SolitarySGN <: InitialData

    η
    v
    label::String
    info::String

    function SolitarySGN(param::NamedTuple; x₀ = 0::Real)
        (η, u, v, mesh) = SolitaryWaveSerreGreenNaghdi(param; x₀)
        init = Init(mesh, η, v)
        label = "Green-Naghdi solitary wave"
        info = "Solitary travelling wave for the Serre-Green-Naghdi model.\n\
        ├─velocity c = $(param.c)\n\
        └─maximum h₀ = $((param.c^2 - 1) / param.ϵ) (from rest state)."

        return new(init.η, init.v, label, info)
    end

    function SolitarySGN(c::Real; ϵ = 1::Real, μ = 1::Real, x₀ = 0::Real, N = 2^12::Int)
        L = 200 / sqrt(3 * (c^2 - 1) / (c^2) / μ) / 2
        xmin, xmax = x₀ - L, x₀ + L
        param = (ϵ = ϵ, μ = μ, c = c, L = L, xmin = xmin, xmax = xmax, N = N)
        sol = SolitarySGN(param; x₀ = x₀)
        return new(sol.η, sol.v, sol.label, sol.info)
    end

end
