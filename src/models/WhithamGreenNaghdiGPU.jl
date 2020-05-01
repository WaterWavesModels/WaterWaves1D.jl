export WhithamGreenNaghdiGPU,mapto,mapfro
using CUDAdrv, CuArrays, CuArrays.CUFFT

"""
    WhithamGreenNaghdiGPU(param;kwargs)

    GPU version of [`WhithamGreenNaghdiGPU`](@ref)

"""
struct WhithamGreenNaghdiGPU <: AbstractModel

    label::String
    datasize::Int
    μ::Real
    ϵ::Real
    x::Vector{Float64}
    F₀::Vector{ComplexF64}
    ∂ₓ::Vector{ComplexF64}
    Π⅔::Vector{Float64}
    Id::BitArray{2}
    FFT::Array{ComplexF64,2}
    IFFTF₀::Array{ComplexF64,2}
    M₀::Array{ComplexF64,2}
    ktol::Real


    function WhithamGreenNaghdiGPU(
        param::NamedTuple;
        SGN = false,
        ktol = 0,
        dealias = 0,
    )

        @show CUDAdrv.name(CuDevice(0))

        if SGN == true
            label = string("Serre-Green-Naghdi")
        else
            label = string("Whitham-Green-Naghdi")
        end
        @info label

        datasize = 2
        μ = param.μ
        ϵ = param.ϵ
        mesh = Mesh(param)
        k = mesh.k
        x = mesh.x
        x₀ = mesh.x[1]

        ∂ₓ = 1im * mesh.k
        F₁ = tanh.(sqrt(μ) * abs.(k)) ./ (sqrt(μ) * abs.(k))
        F₁[1] = 1                 # Differentiation

        if SGN == true
            F₀ = sqrt(μ) * ∂ₓ
        else
            F₀ = 1im * sqrt.(3 * (1 ./ F₁ .- 1)) .* sign.(k)
        end

        K = mesh.kmax * (1 - dealias / (2 + dealias))
        Π⅔ = abs.(mesh.k) .<= K # Dealiasing low-pass filter
        if dealias == 0
            @info "no dealiasing"
            Π⅔ = ones(size(mesh.k))
        else
            @info "dealiasing"
        end

        @info "LU decomposition"

        FFT = exp.(-1im * k * (x .- x₀)')
        IFFT = exp.(1im * k * (x .- x₀)') / length(x)
        IFFTF₀ = IFFT * Diagonal(F₀)
        M₀ = IFFT * Diagonal(F₀) * FFT
        Id = Diagonal(ones(size(x)))

        new(
            label,
            datasize,
            μ,
            ϵ,
            x,
            F₀,
            ∂ₓ,
            Π⅔,
            Id,
            FFT,
            IFFTF₀,
            M₀,
            ktol
        )
    end
end


function (model::WhithamGreenNaghdiGPU)( U::Array{ComplexF64,2} )

    d_fftη = CuArray(U[:, 1])
    d_fftu = CuArray(U[:, 2])
    d_fftv = copy(d_fftu)
    d_h = similar(d_fftu)
    d_u = similar(d_fftu)
    d_hdu = similar(d_fftu)
    d_mo = CuArray(model.M₀)
    d_pi = CuArray(model.Π⅔)
    d_id = CuArray(model.Id)
    d_fo = CuArray(model.F₀)
    d_L = CuArray(model.FFT)
    d_dx = CuArray(model.∂ₓ)
    d_fft = CuArray(model.FFT)
    d_ifft0 = CuArray(model.IFFTF₀)

    fw = CUFFT.plan_fft!(d_fftη)
    bw = CUFFT.plan_ifft!(d_fftv)

    d_fftv .= d_fftu
    d_h .= d_fftη
    bw * d_h
    d_h .*= model.ϵ
    d_h .+= 1

    d_L .= Diagonal(d_h .* d_h .* d_h) * d_ifft0
    d_L .= d_mo * d_L
    d_L .= Diagonal(1 ./ d_h) * d_L
    d_L .= d_fft * d_L
    d_L .= d_id - 1 / 3 * d_L

    d_L, ipiv = CuArrays.CUSOLVER.getrf!(d_L)
    CuArrays.CUSOLVER.getrs!('N', d_L, ipiv, d_fftu)

    d_u .= d_fftu
    d_hdu .= d_u

    bw * d_u

    d_hdu .*= d_fo
    d_hdu .*= d_pi
    bw * d_hdu
    d_hdu .*= d_h
    d_hdu .*= d_hdu
    d_hdu .*= 1 / 2

    bw * d_fftv
    d_fftv .*= d_u
    d_fftv .-= 1 / 2 .* d_u .* d_u
    d_fftv .-= d_hdu

    fw * d_fftv
    d_fftv .*= model.ϵ
    d_fftv .+= d_fftη
    d_fftv .*= d_pi
    d_fftv .*= d_dx

    bw * d_fftη
    d_fftη .*= d_u
    fw * d_fftη
    d_fftη .*= model.ϵ
    d_fftη .+= d_fftu
    d_fftη .*= d_pi
    d_fftη .*= d_dx

    d_fftη .*= -1
    d_fftv .*= -1

    U[:, 1] .= Array(d_fftη)
    U[:, 2] .= Array(d_fftv)

	U[abs.(U) .< model.ktol] .= 0.0

end

"""
    mapto(WhithamGreenNaghdi, data)
`data` is of type `InitialData`, maybe constructed by `Init(...)`.

Performs a discrete Fourier transform with, possibly, dealiasing and Krasny filter.

See documentation of `WhithamGreenNaghdi` for more details.

"""
function mapto(m::WhithamGreenNaghdiGPU, data::InitialData)

    U = [m.Π⅔ .* fft(data.η(m.x)) m.Π⅔ .* fft(data.v(m.x))]
    U[abs.(U).<m.ktol] .= 0
    return U

end

"""
    mapfro(WhithamGreenNaghdi, data)
`data` is of type `Array{Complex{Float64},2}`, e.g. `last(p.data.U)` where `p` is of type `Problem`.

Returns `(η,v,u)`, where
- `η` is the surface deformation;
- `v` is the derivative of the trace of the velocity potential;
- `u` corresponds to the layer-averaged velocity.

Inverse Fourier transform and real part, plus solving the elliptic problem for `u`.

See documentation of `WhithamGreenNaghdi` for more details.
"""
function mapfro(m::WhithamGreenNaghdiGPU, datum::Array{Complex{Float64},2})

    fftη = datum[:, 1]
    h = 1 .+ m.ϵ * ifft(fftη)
    L =
        m.Id - 1 / 3 * m.FFT * Diagonal(1 ./ h) * m.M₀ * Diagonal(h .^ 3) * m.IFFTF₀

    real(ifft(datum[:, 1])), real(ifft(datum[:, 2])), real(ifft(L \ datum[:, 2]))

end
