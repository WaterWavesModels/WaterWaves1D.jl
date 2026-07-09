using WaterWaves1D

export test
"""
    test(param; kwargs...)

Define an object of type `AbstractModel` in view of testing
"""
mutable struct test <: AbstractModel

    label::String
    f!::Function
    g!::Function
    D::AbstractArray
    mapto::Function
    mapfro::Function

    function test(
            param::NamedTuple;
            mesh = Mesh(param),
            label = "test"
        )

        # Set up
        # Pre-allocate useful data
        k = mesh.k
        x = mesh.x

        η = zeros(Complex{Float64}, mesh.N)
        v = zeros(Complex{Float64}, mesh.N)


        # Evolution equations are ∂t U = f(U)
        function f!(U)
            η .= U[1]
            v .= U[2]
            v .*= -1
            #
            U[1] .= v
            return U[2] .= η
            #U .*= 1im
        end
        function g!(U)
            η .= U[1]
            v .= U[2]
            v .*= -1
            #
            U[1] .= v
            return U[2] .= η
        end
        D = [0 * η, 0 * v]
        # Build raw data from physical data.
        # Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
        function mapto(data::InitialData)
            U = [fft(data.η(x)), fft(data.v(x))]
            return U
        end

        # Reconstruct physical variables from raw data
        # Return `(η,v,x)`, where
        # - `η` is the surface deformation;
        # - `v` is the derivative of the trace of the velocity potential;
        # - `x` is the vector of collocation points
        function mapfro(U)
            return real(ifft(U[1])), real(ifft(U[2])), mesh.x
        end

        return new(label, f!, g!, D, mapto, mapfro)
    end
end

"""
    test(param; kwargs...)

Define an object of type `AbstractModel` in view of testing
"""
mutable struct badtest <: AbstractModel

    label::String
    f!::Function
    mapto::Function
    mapfro::Function

    function badtest(
            param::NamedTuple;
            mesh = Mesh(param),
            label = "test"
        )

        # Set up
        # Pre-allocate useful data
        k = mesh.k
        x = mesh.x

        η = zeros(Complex{Float64}, mesh.N)
        v = zeros(Complex{Float64}, mesh.N)


        # Evolution equations are ∂t U = f(U)
        function f!(U)
            a = U[1]
            b = U[2]

            U[1] .= -b
            return U[2] .= a
        end

        # Build raw data from physical data.
        # Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
        function mapto(data::InitialData)
            U = [fft(data.η(x)), fft(data.v(x))]
            return U
        end

        # Reconstruct physical variables from raw data
        # Return `(η,v,x)`, where
        # - `η` is the surface deformation;
        # - `v` is the derivative of the trace of the velocity potential;
        # - `x` is the vector of collocation points
        function mapfro(U)
            return real(ifft(U[1])), real(ifft(U[2])), mesh.x
        end

        return new(label, f!, mapto, mapfro)
    end
end

#---- parameters
param = (
    μ = 0.01,
    ϵ = 1,
    N = 2^6, # number of collocation points
    L = π, # size of the mesh (-L,L)
    T = 1, # final time of computation
    dt = 0.001,  # timestep
    Ns = 1,
);

#---- initial data
ζ(x) = zero(x);
u(x) = -sin.(x);
init = Init(ζ, u);
# ζ et v sont des fonctions.
# Init les mets sous une forme utilisable par le programme


#---- models to compare
using FFTW


using BenchmarkTools

pb = Problem(test(param), init, param, solver = Euler(test(param)), label = "Euler")
pb_naive = Problem(test(param), init, param, solver = Euler_naive(), label = "Euler_naive")
U = pb.data.U[1]

badmodel = badtest(param)
@btime badmodel.f!(U);

model = test(param)
@btime model.f!(U);

solver = EulerExp(test(param))
model = test(param)
pb = Problem(model, init, param, solver = solver, label = "Euler")
solver_naive = EulerExp_naive()
pb_naive = Problem(model, init, param, solver = solver_naive, label = "Euler_naive")

@btime step!(pb.solver, pb.model, U, param.dt);
@btime step!(pb_naive.solver, pb.model, U, param.dt);

@time solve!(pb)
@time solve!(pb_naive)

solver = Euler(test(param))
model = test(param)
pb = Problem(model, init, param, solver = solver, label = "Euler")
solver_naive = Euler_naive()
pb_naive = Problem(test(param), init, param, solver = solver_naive, label = "Euler_naive")

@btime step!(pb.solver, pb.model, U, param.dt);
@btime step!(pb_naive.solver, pb.model, U, param.dt);

@time solve!(pb)
@time solve!(pb_naive)


solver = RK4(test(param))
model = badtest(param)
pb = Problem(model, init, param, solver = solver, label = "RK4")
solverbis = RK4bis(test(param))

pbbis = Problem(model, init, param, solver = solverbis, label = "RK4bis")
solver_naive = RK4_naive()
pb_naive = Problem(model, init, param, solver = solver_naive, label = "Euler_naive")

U = pb.data.U[1]

@btime step!(pb.solver, pb.model, U, param.dt);
@btime step!(pbbis.solver, pb.model, U, param.dt);
@btime step!(pb_naive.solver, pb.model, U, param.dt);


@time solve!(pb)
@time solve!(pbbis)
@time solve!(pb_naive)


using Plots
plot(pb_naive)
