using Test
using FFTW
using LinearAlgebra
using ProgressMeter
using BenchmarkTools
using DeepWaterModels

param = Parameters( ϵ  = 1/2, 
                    N  = 2^12,
                    L  = 10,
                    T  = 5,
                    dt = 0.001)

times   = Times(param.dt, param.T)
solver  = RK4( param.N )
mesh    = Mesh(-param.L, param.L, param.N)
cheng   = Cheng(mesh, param.ϵ)
matsuno = Matsuno(mesh, param.ϵ)
bump    = Bump( param )
problem = Problem( cheng, bump, param, solver )

h  = bump.h
u  = bump.u

@testset "Parameters" begin

    @test param.ϵ  == 0.5
    @test param.N  == 4096
    @test param.L  == 10
    @test param.T  == 5
    @test param.dt == 0.001

end

@testset "Initial data" begin

    @test length(h) == param.N
    @test length(u) == param.N

end

@testset "Problem" begin

    @test true

end

h .= cheng.Pi .* fft(h)
u .= cheng.Pi .* fft(u)

@test !any(isnan,h)
@test !any(isnan,u)

h .= bump.h
u .= bump.u

h .= matsuno.Pi .* fft(h)
u .= matsuno.Pi .* fft(u)

@test !any(isnan,h)
@test !any(isnan,u)
