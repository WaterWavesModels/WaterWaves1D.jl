export bilayerSaintVenant

"""
    bilayerSaintVenant(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for the bilayer
Saint-Venant (or shallow water) model.

# Argument
`param` is of type `NamedTuple` and must contain
- the dimensionless parameter `ϵ` (nonlinearity);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: no dealisasing if set to `0` or `false` (default), standard 3/2 Orszag rule if set to `1` or `true`, otherwise the value sets additionnally a maximal slope of the dealiasing symbol (`2/dealias` models are affected);
- `label`: a label for future references (default is `"Saint-Venant"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `bilayerSaintVenant.f!` to be called in explicit time-integration solvers;
2. a function `bilayerSaintVenant.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `bilayerSaintVenant.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct bilayerSaintVenant <: AbstractModel

    label::String
    f!::Function
    mapto::Function
    mapfro::Function
    info::String

    function bilayerSaintVenant(
            param::NamedTuple;
            mesh = Mesh(param),
            dealias = 0, ktol = 0,
            label = "linearized Saint-Venant",
            γ = 0.9, δ = 1, smooth = false
        )

        # Set up
        ϵ = param.ϵ

        info = "Bilayer Saint-Venant model.\n"
        info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
        info *= "├─density ratio parameter γ=$(param.γ).\n"
        info *= "├─depth ratio parameter δ=$(param.δ).\n"
        if dealias == 0
            info *= "└─No dealiasing. "
        else
            info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity: \n"
            if smooth
                info *= "└─Lipschitz cut-off. "
            else
                info *= "└─Sharp cut-off. "
            end
        end
        if ktol == 0
            info *= "No Krasny filter. "
        else
            info *= "Krasny filter with tolerance $ktol."
        end
        info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

        # Pre-allocate useful data
        k = mesh.k
        x = mesh.x
        dk = mesh.dk
        ∂ₓ = 1im * k

        if dealias == 0
            Π⅔ = ones(size(k)) # no dealiasing (Π⅔=Id)
            Π = ones(size(k)) # no dealiasing (Π⅔=Id)
        else
            K = (mesh.kmax - mesh.kmin) / (3 * dealias)
            Π⅔ = abs.(k) .<= K # Dealiasing low-pass filter
            if smooth == 0
                Π = abs.(k) .<= K # Dealiasing low-pass filter
            else
                Π = max.(0, min.(1, 2 / smooth * (1 .- abs.(k) / K))) .^ 2
            end
        end
        η = zeros(Float64, mesh.N)
        v = zeros(Float64, mesh.N)
        h1 = zeros(Float64, mesh.N)
        h2 = zeros(Float64, mesh.N)


        # Evolution equations are ∂t U = f(U)
        function f!(U)
            #fftη .= U[1]
            #fftv .= U[2]
            η .= real(ifft(U[1]))
            h1 .= 1 .- ϵ * η
            h2 .= 1 / δ .+ ϵ * η
            v .= real(ifft(U[2]))

            U[1] .= -∂ₓ .* Π .* fft(h1 .* h2 .* v ./ (h1 .+ γ * h2))
            U[2] .= -∂ₓ .* Π .* fft(η .+ ϵ / 2 * ((h1 .^ 2 - γ * h2 .^ 2) .* v .^ 2) ./ ((h1 .+ γ * h2) .^ 2))
            for u in U
                u[abs.(u) .< ktol] .= 0
            end
            return
        end

        # Build raw data from physical data.
        # Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
        function mapto(data::InitialData)
            U = [Π⅔ .* fft(data.η(x)), Π⅔ .* fft(data.v(x))]
            for u in U
                u[abs.(u) .< ktol] .= 0
            end
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

        return new(label, f!, mapto, mapfro, info)
    end
end
