using FFTW
using WaterWaves1D
using RecursiveArrayTools
using TimerOutputs

const to = TimerOutput()

struct TestSaintVenant2D <: AbstractModel

    label::String
    f!::Function
    mapto::Function
    mapfro::Function

    function TestSaintVenant2D(
        param::NamedTuple;
        mesh = Mesh(param),
        dealias = 1,
        ktol = 0,
        smooth = false,
        hamiltonian = false,
        label = "Saint-Venant",
    )

        ϵ = param.ϵ

        kx = mesh.k
        ky = mesh.k
        x = mesh.x
        y = mesh.x
        nx = mesh.N
        ny = mesh.N

        K = (mesh.kmax - mesh.kmin) / 3
        Πx⅔ = abs.(kx) .<= K
        Πy⅔ = abs.(ky') .<= K
        Πx = abs.(kx) .<= K
        Πy = abs.(ky') .<= K


        η = zeros(ComplexF64, (nx, ny))
        vx = zeros(ComplexF64, (nx, ny))
        vy = zeros(ComplexF64, (nx, ny))
        fftη = zeros(ComplexF64, (nx, ny))
        fftvx = zeros(ComplexF64, (nx, ny))
        fftvy = zeros(ComplexF64, (nx, ny))

        Ix = similar(fftη)
        Iy = similar(fftη)
        I = similar(fftη)
        Jx = similar(fftη)
        Jy = similar(fftη)

        ϵΠx = (√ϵ) * Πx
        ϵΠy = (√ϵ) * Πy
        ∂y = -1im * ky'
        ∂x = -1im * kx

        function f!(U)
            fftη .= U[1]
            fftvx .= U[2]
            fftvy .= U[3]

            η .= ifft(fftη)
            vx .= ifft(fftvx)
            vy .= ifft(fftvy)

            Jx .= fftvx
            Jx .*= ∂x
            my_ifft!(Ix, Jx, I)
            Ix .*= vx

            Jy .= fftvx
            Jy .*= ∂y
            my_ifft!(Iy, Jy, I)
            Iy .*= vy

            Ix .+= Iy
            my_fft!(I, Ix, Iy)

            I .*= ϵΠx
            I .*= ϵΠy

            U[2] .= fftη
            U[2] .*= ∂x
            U[2] .+= I

            Jx .= fftvy
            Jx .*= ∂x
            my_ifft!(Ix, Jx, I)
            Ix .*= vx

            Jy .= fftvy
            Jy .*= ∂y
            my_ifft!(Iy, Jy, I)
            Iy .*= vy

            Ix .+= Iy
            my_fft!(I, Ix, Iy)

            I .*= ϵΠx
            I .*= ϵΠy

            U[3] .= fftη
            U[3] .*= ∂y
            U[3] .+= I

            vx .*= η
            Ix .= fft(vx)
            Ix .*= ϵΠx
            Ix .*= ϵΠy
            fftvx .+= Ix
            fftvx .*= ∂x


            vy .*= η
            Iy .= fft(vy)
            Iy .*= ϵΠy
            Iy .*= ϵΠx
            fftvy .+= Iy
            fftvy .*= ∂y

            fftvx .+= fftvy
            U[1] .= fftvx

        end

        function mapto(data::InitialData)
            U = [
                Πx⅔ .* Πy⅔ .* fft(data.η(x, y)),
                Πx⅔ .* Πy⅔ .* fft(data.vx(x, y)),
                Πx⅔ .* Πy⅔ .* fft(data.vy(x, y)),
            ]

            for u in U
                u[abs.(u).<ktol] .= 0
            end
            return U
        end

        function mapfro(U)
            real(ifft(U[1])), real(ifft(U[2])), real(ifft(U[3])), x, y
        end

        new("Saint-Venant", f!, mapto, mapfro)
    end
end

param = ( μ = 0.01, ϵ = 0.1, N = 2^8, L = π, T = 0.1, dt = 0.001, Ns = 1)

K = floor(param.N / 2 * 2 / 3) - 1

ζ(x, y) = 1 / 2 * cos.(x) .* cos.(y)'
ux(x, y) = cos.(y') .* sin.(x) + cos.(K * y') .* sin.(K * x) / K
uy(x, y) = -sin.(y') .* cos.(x) - sin.(K * y') .* cos.(K * x) / K

init = Init2D(ζ, ux, uy)

model1 = SaintVenant2D(param; dealias = 1)
model2 = SaintVenant2D_fast( param; dealias = 1)
model3 = SaintVenant2D_fast( param; dealias = 1, large_data = true)
model4 = TestSaintVenant2D(param)


solver = RK4(model1.mapto(init))


problem1 = Problem(model1, init, param; solver = solver)
problem2 = Problem(model2, init, param; solver = solver)
problem3 = Problem(model3, init, param; solver = solver)

#---- Solve problems

@timeit to "classic problem" solve!(problem1)
@timeit to "fast problem" solve!(problem2)
@timeit to "fast problem with large data" solve!(problem3)


show(to)
