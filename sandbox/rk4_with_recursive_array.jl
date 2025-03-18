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
problem2 = Problem(model2, init, param; solver = solver)
problem3 = Problem(model3, init, param; solver = solver)
problem4 = Problem(model4, init, param; solver = solver)

#---- Solve problems

@timeit to "classic problem" solve!(problem1)
@timeit to "fast problem" solve!(problem2)
@timeit to "fast problem with large data" solve!(problem3)
@timeit to "test problem " solve!(problem4)


# Sanity check
η1, vx1, vy1, x1, y1 = problem1.model.mapfro(problem1.data.U[end])
η2, vx2, vy2, x2, y2 = problem2.model.mapfro(problem2.data.U[end])
η3, vx3, vy3, x3, y3 = problem3.model.mapfro(problem3.data.U[end])
η4, vx4, vy4, x4, y4 = problem4.model.mapfro(problem4.data.U[end])


@show norm(η1 - η2), norm(vx1 - vx2), norm(vy1 - vy2)
@show norm(η1 - η3), norm(vx1 - vx3), norm(vy1 - vy3)
@show norm(η1 - η4), norm(vx1 - vx4), norm(vy1 - vy4)

show(to)

#=

problem_fast ≈ problem_fast2,problem_fast ≈ problem
nothing
# # Profiling allocations to improve the code
# using BenchmarkTools
# using Profile,PProf

# U=problem_fast.data.U[1]
# step!(solver, modelfast, U, 1)

# Profile.Allocs.clear()
# Profile.Allocs.@profile sample_rate=1 step!(solver, modelfast, U, 1)
# PProf.Allocs.pprof(from_c=false)



#---- Solve other problems



model_hamiltonian=SaintVenant2D_fast(param; dealias = 1, hamiltonian = true)
model_smooth=SaintVenant2D_fast(param; dealias = 1, smooth = true, hamiltonian = false)

problem_hamiltonian = Problem(  model_hamiltonian, init, param ; solver=solver )
problem_smooth = Problem(  model_smooth, init, param ; solver=solver )

@time solve!(problem_hamiltonian)
@time solve!(problem_smooth)

#---- Visualization

η,vx,vy, = solution(problem_fast;T=1)
hη,hvx,hvy, = solution(problem_hamiltonian;T=1)
sη,svx,svy, = solution(problem_smooth;T=1)


mesh=Mesh(param)
kx=mesh.k;ky=mesh.k;x=mesh.x;y=mesh.x;
∂y = 1im * ky' ;
∂x = 1im .* kx;

function ∂fy( f )
    real.(ifft(∂y .* fft(f, 2), 2))
end
function ∂fx( f )
    real.(ifft(∂x .* fft(f, 1), 1))
end


plot(x,y,vx,st=:surface)
plot(x,y,hvx,st=:surface)
plot(x,y,svx,st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(vx)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(hvx)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(svx)))),st=:surface)

plot(x,y,η-0*sη,st=:surface)
plot(x,y,hη,st=:surface)
plot(x,y,sη,st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(η)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(hη)))),st=:surface)
plot(fftshift(kx),fftshift(ky),log10.(abs.(fftshift(ifft(sη)))),st=:surface)

plot(x,y,∂fx(vx),st=:surface)

plot(x,y,∂fx(∂fx(vx)),st=:surface)
plot(x,y,∂fx(∂fx(hvx)),st=:surface)
plot(x,y,∂fx(∂fx(svx)),st=:surface)


@gif for time in LinRange(0,1,101)
    η,vx,vy, = solution(problem;T=time)
    plot(x,y,η,st=:surface,zlims=(-1,2));
end
# anim = @animate for time in LinRange(0,param.T,101)
#     η,vx,vy, = sol(problem;T=time)
#     plot(x,y,η,st=:surface);
#     zlims!(-1, 2)
# end
# gif(anim,"SaintVenant2D.gif")

# Energy preservation
function energy(η,vx,vy)
    (η.^2 .+ (1 .+ η) .* (vx.^2 .+ vy.^2))/2
end

function energy_norm(η,vx,vy)
    real.(fft(energy(η,vx,vy)))[1,1]/length(η)
end

function energy_norm(problem,time)
    η,vx,vy, = solution(problem;T=time)
    energy_norm(η,vx,vy)
end

times=LinRange(0,1,11)
plot(times,[energy_norm(problem,t) for t in times])
plot!(times,[energy_norm(problem_hamiltonian,t) for t in times])
plot!(times,[energy_norm(problem_smooth,t) for t in times])



=#
