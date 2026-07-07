export CnoidalWaveSerreGreenNaghdi, CnoidalSGN
using Elliptic

"""
    CnoidalWaveSerreGreenNaghdi(param; P=1)

Compute the Serre-Green-Naghdi cnoidal wave with prescribed `h₀<h₁<h₂`.
`h₁` is the minimum, `h₂` is the maximum of the wave.
As `h₀ -> h₁`, the cnoidal wave converges towards the solitary wave.
See for instance [GavrilyukNkongaShyueEtAl2020](@citet).

# Arguments
- `param :: NamedTuple`: parameters of the problem containing `h₀<h₁<h₂` and dimensionless parameters `ϵ` and `μ`, and number of collocation points `N`.
- `P :: Int`: (keyword, optional, default = 1) the number of periods of the cnoidal wave in the constructed mesh.

# Return values
`(η,u,v,mesh,param)` with
- `η :: Vector{Float64}`: surface deformation;
- `u :: Vector{Float64}`: layer-averaged velocity;
- `v :: Vector{Float64}`: derivative of the trace of the velocity potential at the surface;
- `mesh :: Mesh`: mesh collocation points;
- `param :: NamedTuple`: useful parameters
"""
function CnoidalWaveSerreGreenNaghdi(
        param::NamedTuple;
        P = 1::Int
    )

    ϵ = param.ϵ
    μ = param.μ

    h₀ = param.h₀
    h₁ = param.h₁
    h₂ = param.h₂
    c = sqrt(h₀ * h₁ * h₂)
    m = sqrt((h₂ - h₁) / (h₂ - h₀))
    κ = sqrt(3 * (h₂ - h₀)) / (2 * c) / sqrt(μ)
    λ = Elliptic.K(m^2) / κ
    mesh = Mesh((L = P * λ, N = param.N))
    formula = h₁ .- 1 .+ (h₂ - h₁) * (Jacobi.cn.(κ * mesh.x, m^2) .^ 2)

    a₀ = h₀
    a₁ = h₂ - h₀
    formula2 = a₀ .- 1 .+ a₁ * (Jacobi.dn.(κ * mesh.x, m^2) .^ 2)
    H₀ = a₀ + a₁ * Elliptic.E(m^2) / Elliptic.K(m^2)
    u2 = c * (1 ./ H₀ .- 1 ./ (1 .+ formula2))
    param = (h₀ = h₀, h₁ = h₁, h₂ = h₂, a₀ = a₀, a₁ = a₁, H₀ = H₀, c = c, λ = λ, m = m, κ = κ)


    η = formula / ϵ
    h = 1 .+ ϵ * η
    u = c * η ./ h
    k = mesh.k
    Dx = 1im * k
    F₀ = sqrt(μ) * Dx
    DxF(v) = real.(ifft(F₀ .* fft(v)))
    v = u - 1 / 3 ./ h .* (DxF(h .^ 3 .* DxF(u)))
    #
    # h2 = 1 .+ ϵ*formula2
    # v2 = u2 - 1/3 ./h2 .* (DxF(h2.^3 .*DxF(u2)))

    return (η, u, v, mesh, param)

end

"""
    CnoidalSGN(param; P=1)

Build the initial data associated with `CnoidalWaveSerreGreenNaghdi(param; P=1)`, of type `InitialData`,
to be used in initial-value problems `Problem(model, initial::InitialData, param)`.
"""
struct CnoidalSGN <: InitialData

    η
    v
    label::String
    info::String

    function CnoidalSGN(param; P = 1)
        (η, u, v, mesh, para) = CnoidalWaveSerreGreenNaghdi(param; P)
        init = Init(mesh, η, v)
        label = "Green-Naghdi cnoidal wave"
        info = "Cnoidal travelling wave for the Serre-Green-Naghdi model.\n\
        ├─velocity c = $(para.c)\n├───period P = $(2(para.λ))\n\
        ├─maximum h₂ = $(para.h₂) (from bottom)\n└─minimum h₁ = $(para.h₁) (from bottom)."

        return new(init.η, init.v, label, info)
    end
end
