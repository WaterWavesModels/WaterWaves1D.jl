using FFTW
using LinearAlgebra
using TimerOutputs
using WaterWaves1D

const to = TimerOutput()

param = (μ = 0.01, ϵ = 0.1, N = 2^8, L = π, T = 0.1, dt = 0.001, Ns = 1)

K = floor(param.N / 2 * 2 / 3) - 1

ζ(x, y) = 1 / 2 * cos.(x) .* cos.(y)'
ux(x, y) = cos.(y') .* sin.(x) + cos.(K * y') .* sin.(K * x) / K
uy(x, y) = -sin.(y') .* cos.(x) - sin.(K * y') .* cos.(K * x) / K

init = Init2D(ζ, ux, uy)

model1 = SaintVenant2D(param; dealias = 1)
model2 = SaintVenant2D_fast(param; dealias = 1)
model3 = SaintVenant2D_fast(param; dealias = 1, large_data = true)
model4 = TestSaintVenant2D(param)


solver = RK4(model1.mapto(init))


problem1 = Problem(model1, init, param; solver = solver)
# problem2 = Problem(model2, init, param; solver = solver)
# problem3 = Problem(model3, init, param; solver = solver)

solver = TestRK4(model4.mapto(init))

problem4 = Problem(model4, init, param; solver = solver)

#---- Solve problems

@timeit to "classic problem" solve!(problem1)
# @timeit to "fast problem" solve!(problem2)
# @timeit to "fast problem with large data" solve!(problem3)
@timeit to "test problem " solve!(problem4)


# Sanity check
η1, vx1, vy1, x1, y1 = problem1.model.mapfro(problem1.data.U[end])
# η2, vx2, vy2, x2, y2 = problem2.model.mapfro(problem2.data.U[end])
# η3, vx3, vy3, x3, y3 = problem3.model.mapfro(problem3.data.U[end])
η4, vx4, vy4, x4, y4 = problem4.model.mapfro(problem4.data.U[end])


# @show norm(η1 - η2), norm(vx1 - vx2), norm(vy1 - vy2)
# @show norm(η1 - η3), norm(vx1 - vx3), norm(vy1 - vy3)
@show norm(η1 - η4), norm(vx1 - vx4), norm(vy1 - vy4)

show(to)

