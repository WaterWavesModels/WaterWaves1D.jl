using Pkg

GPU_ENABLED = haskey(Pkg.installed(), "CUDAdrv")

if GPU_ENABLED

    using CUDAdrv, CuArrays, CuArrays.CUFFT

end


export WhithamGreenNaghdiGPU

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
    F₀::Vector{Complex{Float64}}
    ∂ₓ::Vector{Complex{Float64}}
    Π⅔::Vector{Float64}
    Id::BitArray{2}
    FFT::Array{Complex{Float64},2}
    IFFT::Array{Complex{Float64},2}
    IFFTF₀::Array{Complex{Float64},2}
    M₀::Array{Complex{Float64},2}
    h::Vector{Complex{Float64}}
    u::Vector{Complex{Float64}}
    fftv::Vector{Complex{Float64}}
    fftη::Vector{Complex{Float64}}
    fftu::Vector{Complex{Float64}}
    hdu::Vector{Complex{Float64}}
    L::Array{Complex{Float64},2}


    function WhithamGreenNaghdiGPU(
        param::NamedTuple;
        iterate = true,
        SGN = false,
        dealias = 0,
    )

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
        if iterate == true
            @info "GMRES method"
        else
            @info "LU decomposition"
        end
        FFT = exp.(-1im * k * (x .- x₀)')
        IFFT = exp.(1im * k * (x .- x₀)') / length(x)
        M₀ = IFFT * Diagonal(F₀) * FFT
        IFFTF₀ = IFFT * Diagonal(F₀)
        Id = Diagonal(ones(size(x)))
        h = zeros(Complex{Float64}, mesh.N)
        u, fftv, fftη, fftu, hdu = (similar(h),) .* ones(5)
        L = similar(FFT)

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
            IFFT,
            IFFTF₀,
            M₀,
            h,
            u,
            fftv,
            fftη,
            fftu,
            hdu,
            L,
        )
    end
end


function (m::WhithamGreenNaghdi)( U, dt )

    d_fftη = CuArray(U[:, 1])
    d_fftu = CuArray(U[:, 2])
    d_fftv = copy(d_fftu)
    d_h = CuArray(model.h)
    d_u = CuArray(model.u)
    d_hdu = CuArray(model.hdu)
    d_fft = CuArray(model.FFT)
    d_ifft = CuArray(model.IFFTF₀)
    d_mo = CuArray(model.M₀)
    d_pi = CuArray(model.Π⅔)
    d_id = CuArray(model.Id)
    d_fo = CuArray(model.F₀)
    d_L = CuArray(model.L)
    d_dx = CuArray(model.∂ₓ)

    d_fftv .= d_fftu

    fw = CUFFT.plan_fft!(d_fftη)
    bw = CUFFT.plan_ifft!(d_fftv)

    d_h .= d_fftη
    bw * d_h
    d_h .*= m.ϵ
    d_h .+= 1

    d_L .= Diagonal(d_h .* d_h .* d_h) * d_ifft
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
    d_fftv .*= m.ϵ
    d_fftv .+= d_fftη
    d_fftv .*= d_pi
    d_fftv .*= d_dx

    bw * d_fftη
    d_fftη .*= d_u
    fw * d_fftη
    d_fftη .*= m.ϵ
    d_fftη .+= d_fftu
    d_fftη .*= d_pi
    d_fftη .*= d_dx

    d_fftη .*= -1
    d_fftv .*= -1

    U[:, 1] .= Array(d_fftη)
    U[:, 2] .= Array(d_fftv)

    U[abs.(U).<m.ktol] .= 0.0

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

    m.fftη .= datum[:, 1]
    m.h .= 1 .+ m.ϵ * ifft(m.fftη)
    m.L .=
        m.Id - 1 / 3 * m.FFT * Diagonal(1 ./ m.h) * m.M₀ * Diagonal(m.h .^ 3) * m.IFFTF₀

    real(ifft(datum[:, 1])), real(ifft(datum[:, 2])), real(ifft(m.L \ datum[:, 2]))

end
