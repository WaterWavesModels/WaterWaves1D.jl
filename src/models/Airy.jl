export Airy
@doc raw"""
	Airy(param; mesh, label)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the linear ([Airy](https://en.wikipedia.org/wiki/Airy_wave_theory)) water waves equations:
```math
  \left\{\begin{array}{l}
  ∂_tη-\tfrac{1}{\sqrtμ ν} T^μv =0,\\[1ex]
  ∂_tv+∂_xη=0,
  \end{array}\right.
```
where ``η`` is the surface deformation, ``v=∂_xψ`` is the derivative of the trace of the velocity potential at the surface, and
```math
T^μ=-{\rm i}\tanh(\sqrtμ D)
```
is the Fourier multiplier sometimes called "Tilbert transform"
(related to the [Hilbert transform](https://en.wikipedia.org/wiki/Hilbert_transform#Relationship_with_the_Fourier_transform), the latter arising in the infinite layer configuration, ``μ=∞``).

# Arguments
- `param::NamedTuple` must contain
    - the shallowness parameter `μ` (set the infinite-layer case if `μ=Inf`);
    - optionally, `ν` the shallow/deep water scaling factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise;
    - numerical parameters to construct the mesh of collocation points, if `mesh` is not provided.
- `mesh  :: Mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`.
- `label :: String`: a label for future references (default is `"linear (Airy)"`).


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `Airy.f!` to be called in explicit time-integration solvers;
2. a function `Airy.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `Airy.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct Airy <: AbstractModel

    label::String
    f!::Function
    mapto::Function
    mapfro::Function
    info::String


    function Airy(
            param::NamedTuple; # param is a NamedTuple containing all necessary parameters
            mesh = Mesh(param),
            IL = false,
            label = "linear (Airy)"  # using a keyword argument allows the user to supersede the default label.
        )

        # Set up
        μ = param.μ
        if !in(:ν, keys(param))
            if μ > 1
                ν = 1 / sqrt(μ)
                nu = "1/√μ (deep water case)"
            else
                ν = 1
                nu = "1 (shallow water case)"
            end
        else
            ν = param.ν
            nu = "$ν"
        end
        if μ == Inf || ν == 0 || IL == true # infinite layer case
            IL = true   # IL (=Infinite layer) is a flag to be used thereafter
            μ = 1; ν = 1  # Then we should set μ=ν=1 in subsequent formula.
        end

        # Print information
        info = "$label model.\n"
        if IL == true
            info *= "└─Infinite depth case.\n"
        else
            info *= "└─Shallowness parameter μ=$μ, \
                    scaling parameter ν=$nu.\n"
        end
        info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

        # Collocation points and Fourier modes
        x, k = mesh.x, mesh.k
        # Fourier multipliers
        ∂ₓ = 1im * k            # Differentiation
        if IL == true
            ∂ₓF₁ = 1im * sign.(k)
        else
            ∂ₓF₁ = 1im * tanh.(√μ * k) / (√μ * ν)
        end
        # Pre-allocation
        fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)

        # Evolution equations are ∂t U = f(U)
        function f!(U)

            fftη .= U[1]
            fftv .= U[2]

            U[1] .= -∂ₓF₁ .* fftv
            return U[2] .= -∂ₓ .* fftη

        end

        # Build raw data from physical data (discrete Fourier transform)
        function mapto(data::InitialData)
            return [fft(data.η(x)), fft(data.v(x))]
        end

        # Reconstruct physical variables from raw data
        # Return `(η,v,x)`, where
        # - `η` is the surface deformation;
        # - `v` is the derivative of the trace of the velocity potential;
        # - `x` is the vector of collocation points
        function mapfro(U)
            return real(ifft(U[1])), real(ifft(U[2])), mesh.x
        end

        return new(label, f!, mapto, mapfro, info)

    end
end
