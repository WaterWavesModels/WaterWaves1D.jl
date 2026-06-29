export SaintVenant2D, SaintVenant2D_fast

@doc raw"""
    SaintVenant2D(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for the
[Saint-Venant1871](@citet) (or shallow water) model
```math
  \left\{\begin{array}{l}
  ∂_tη+\nabla\cdot((1+ϵη)\bm{v})=0,\\[1ex]
  ∂_tv+\nabla η+ϵ(\bm{v}\cdot\nabla)\bm{v}=0,
  \end{array}\right.
```
where ``η`` is the surface deformation and ``v`` the horizontal velocity (essentially constant along the vertical variable).

# Argument
`param` is of type `NamedTuple` and must contain
- the dimensionless parameter `ϵ` (nonlinearity);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: no dealisasing if set to `0` or `false` (default), otherwise `1/(3*dealias)` modes are set to `0` (corresponding to standard 2/3 Orszag rule if `dealias` is set to `1` or `true`);
- `smooth`: A smooth low-pass filter (whose scaling is defined by ) if set to `0` or `false` (default), otherwise only `2/(3*dealias)*(1-smooth/2)` modes are kept untouched;
- `label`: a label for future references (default is `"Saint-Venant"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `SaintVenant2D.f!` to be called in explicit time-integration solvers;
2. a function `SaintVenant2D.mapto` which from `(η,vx,vy)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `SaintVenant2D.mapfro` which from such data matrix returns the Tuple of real vectors `(η,vx,vy,x,y)`, where
    - `η` is the values of surface deformation at collocation points `(x,y)`;
    - `vx,vy` are the velocity fields at collocation points `(x,y)`.

See also [`SaintVenant2D_fast`](@ref).

"""
mutable struct SaintVenant2D <: AbstractModel

    label::String
    f!::Function
    mapto::Function
    mapfro::Function
    info::String

    function SaintVenant2D(
        param::NamedTuple;
        mesh = Mesh(param),
        dealias = 0,
        ktol = 0,
        smooth = false,
        hamiltonian = false,
        label = "Saint-Venant",
    )

        # Set up
        ϵ = param.ϵ

        if hamiltonian
            info = "2D Saint-Venant model (hamiltonian).\n"
        else
            info = "2D Saint-Venant model.\n"
        end
        info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
        if dealias == 0
            info *= "└─No dealiasing. "
        else
            info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity: \n"
            if dealias == 1
                info *= "└─Sharp cut-off. "
            else
                info *= "└─Lipschitz cut-off. "
            end
        end
        if ktol == 0
            info *= "No Krasny filter. "
        else
            info *= "Krasny filter with tolerance $ktol."
        end
        info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

        # Pre-allocate useful data
        kx = mesh.k
        ky = mesh.k
        x = mesh.x
        y = mesh.x
        nx = mesh.N
        ny = mesh.N
        ∂y = 1im * ky'
        ∂x = 1im .* kx

        if dealias == 0
            Π⅔ = ones(nx, ny) # no dealiasing (Π⅔=Id)
            Πx = ones(nx, 1) # no dealiasing (Π⅔=Id)
            Πy = ones(1, ny) # no dealiasing (Π⅔=Id)
        else
            K = (mesh.kmax - mesh.kmin) / (3 * dealias)
            Π⅔ = (abs.(ky') .<= K) .* (abs.(kx) .<= K) # Dealiasing low-pass filter
            if smooth == 0
                Πx = abs.(kx) .<= K # Dealiasing low-pass filter
                Πy = abs.(ky') .<= K # Dealiasing low-pass filter
            else
                Πx = max.(0, min.(1, 2 / smooth * (1 .- abs.(kx) / K))) .^ 2
                Πy = max.(0, min.(1, 2 / smooth * (1 .- abs.(ky') / K))) .^ 2
            end
            #P   = min.(( abs.(mesh.k) .<= K).*( (-abs.(mesh.k) .+ K)*dealias/mesh.dk ),1)	
            #P  =min.(( abs.(mesh.k) .<= K).*( (-abs.(mesh.k) .+ K)*dealias/mesh.dk .+1/2 ),1)	
        end

        η = zeros(Float64, (nx, ny))
        vx = zeros(Float64, (nx, ny))
        vy = zeros(Float64, (nx, ny))
        fftη = zeros(Complex{Float64}, (nx, ny))
        fftvx = zeros(Complex{Float64}, (nx, ny))
        fftvy = zeros(Complex{Float64}, (nx, ny))

        # Evolution equations are ∂t U = f(U)
        function f!(U)
            fftη .= U[1]
            fftvx .= U[2]
            fftvy .= U[3]

            η .= real(ifft(fftη))
            vx .= real(ifft(fftvx))
            vy .= real(ifft(fftvy))

            U[1] .=
                -∂x .* (fftvx .+ ϵ * Πx .* Πy .* fft(η .* vx)) +
                -∂y .* (fftvy .+ ϵ * Πx .* Πy .* fft(η .* vy))
            if hamiltonian
                U[2] .= -∂x .* (fftη .+ ϵ / 2 * Πx .* Πy .* fft(vx .^ 2 + vy .^ 2))
                U[3] .= -∂y .* (fftη .+ ϵ / 2 * Πx .* Πy .* fft(vx .^ 2 + vy .^ 2))
            else
                U[2] .=
                    -∂x .* (fftη) .-
                    ϵ * Πx .* Πy .* fft(vx .* ifft(∂x .* fftvx) + vy .* ifft(∂y .* fftvx))
                U[3] .=
                    -∂y .* (fftη) .-
                    ϵ * Πx .* Πy .* fft(vx .* ifft(∂x .* fftvy) + vy .* ifft(∂y .* fftvy))
            end
            for u in U
                u[abs.(u).<ktol] .= 0
            end
        end

        # Build raw data from physical data.
        # Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
        function mapto(data::InitialData)
            #U = cat(Π⅔.* fft(data.η(x,y)), Π⅔.*fft(data.vx(x,y)), Π⅔.*fft(data.vy(x,y)),dims=3)
            U = [
                Π⅔ .* fft(data.η(x, y)),
                Π⅔ .* fft(data.vx(x, y)),
                Π⅔ .* fft(data.vy(x, y)),
            ]

            for u in U
                u[abs.(u).<ktol] .= 0
            end
            return U
        end

        # Reconstruct physical variables from raw data
        # Return `(η,vx,vy,x,y)`, where
        # - `η` is the surface deformation;
        # - `vx` is the horizontal velocity field;
        # - `vy` is the horizontal velocity field;
        # - `(x,y)` is the matrix of collocation points
        function mapfro(U)
            real(ifft(U[1])), real(ifft(U[2])), real(ifft(U[3])), x, y
        end

        new(label, f!, mapto, mapfro, info)
    end
end

"""
	SaintVenant2D_fast(param;kwargs)

Same as [`SaintVenant2D`](@ref), but faster.

If the optional argument `large_data` is set to `true` (default is `false`), 
then the standard `fft` and `ifft` functions (instead of plans) are used.
This may be faster, while using more allocations.
"""
mutable struct SaintVenant2D_fast <: AbstractModel

    label::String
    f!::Function
    mapto::Function
    mapfro::Function
    info::String

    function SaintVenant2D_fast(
        param::NamedTuple;
        mesh = Mesh(param),
        dealias = 0,
        ktol = 0,
        smooth = false,
        hamiltonian = false,
        large_data = false,
        label = "Saint-Venant",
    )

        # Set up
        ϵ = param.ϵ

        if hamiltonian
            info = "2D Saint-Venant model (hamiltonian).\n"
        else
            info = "2D Saint-Venant model.\n"
        end
        info *= "├─Nonlinearity parameter ϵ=$(param.ϵ).\n"
        if dealias == 0
            info *= "└─No dealiasing. "
        elseif dealias == 1
            info *= "├─Dealiasing with Orszag's rule adapted to power 2 nonlinearity. \n"
        else
            info *= "├─Dealiasing at user-defined mode. \n"
        end
        if smooth != 0
            info *= "└─Lipschitz low-pass filter is applied. "
        end
        if ktol == 0
            info *= "No Krasny filter. "
        else
            @error("The Krasny filter will *not* be applied. Use 'SaintVenant' instead")
        end
        info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

        # Pre-allocate useful data
        kx = mesh.k
        ky = mesh.k
        x = mesh.x
        y = mesh.x
        nx = mesh.N
        ny = mesh.N

        if dealias == 0 # no dealiasing (Π⅔=Id)
            Πx⅔ = ones(nx, 1)
            Πy⅔ = ones(1, ny)
            Πx = ones(nx, 1)
            Πy = ones(1, ny)

        else
            # Dealiasing low-pass filter
            K = (mesh.kmax - mesh.kmin) / (3 * dealias)
            Πx⅔ = abs.(kx) .<= K
            Πy⅔ = abs.(ky') .<= K
            # Smooth or sharp low-pass filter
            if smooth == 0
                Πx = abs.(kx) .<= K
                Πy = abs.(ky') .<= K
            else
                Πx = max.(0, min.(1, 2 / smooth * (1 .- abs.(kx) / K))) .^ 2
                Πy = max.(0, min.(1, 2 / smooth * (1 .- abs.(ky') / K))) .^ 2
            end
        end


        η = zeros(Complex{Float64}, (nx, ny))
        vx = zeros(Complex{Float64}, (nx, ny))
        vy = zeros(Complex{Float64}, (nx, ny))
        fftη = zeros(Complex{Float64}, (nx, ny))
        fftvx = zeros(Complex{Float64}, (nx, ny))
        fftvy = zeros(Complex{Float64}, (nx, ny))
        Ix = similar(fftη)
        Iy = similar(fftη)
        I = similar(fftη)
        Jx = similar(fftη)
        Jy = similar(fftη)

        ϵΠx = (√ϵ) * Πx
        ϵΠy = (√ϵ) * Πy
        ∂y = -1im * ky'
        ∂x = -1im * kx

        FFTW.set_num_threads(4)
        Px = plan_fft(η, 1)#, flags=FFTW.PATIENT)    
        #Py = plan_fft(fᵗ, 1, flags=FFTW.PATIENT)
        Py = plan_fft(η, 2)#, flags=FFTW.PATIENT)    
        iPx = plan_ifft(η, 1)#, flags=FFTW.PATIENT)    
        #Py = plan_ifft(fᵗ, 1, flags=FFTW.PATIENT)
        iPy = plan_ifft(η, 2)#, flags=FFTW.PATIENT)    


        function my_ifft!(f, f̂, cache; large_data = large_data)

            if large_data
                f .= ifft(f̂)
            else
                mul!(cache, iPx, f̂)
                mul!(f, iPy, cache)
                #ldiv!(cache, Px, f̂ )
                #ldiv!(f, Py, cache )

            end
        end
        function my_fft!(f̂, f, cache; large_data = large_data)
            if large_data
                f̂ .= fft(f)
            else
                mul!(cache, Px, f)
                mul!(f̂, Py, cache)

            end
        end

        # Evolution equations are ∂t U = f(U)
        function f!(U)
            fftη .= U[1]
            fftvx .= U[2]
            fftvy .= U[3]

            my_ifft!(η, fftη, I)
            my_ifft!(vx, fftvx, I)
            my_ifft!(vy, fftvy, I)

            if hamiltonian
                Jx .= vx
                Jx .^= 2
                Jy .= vy
                Jy .^= 2
                Jx .+= Jy
                Jx ./= 2
                my_fft!(I, Jx, Ix)  # I = 1/2(vx^2+vy^2) in Fourier space

                I .*= ϵΠx
                I .*= ϵΠy
                I .+= fftη

                U[2] .= I
                U[2] .*= ∂x

                U[3] .= I
                U[3] .*= ∂y

            else
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
            end

            vx .*= η
            my_fft!(Ix, vx, I)
            Ix .*= ϵΠx
            Ix .*= ϵΠy
            fftvx .+= Ix
            fftvx .*= ∂x


            vy .*= η
            my_fft!(Iy, vy, I)
            Iy .*= ϵΠy
            Iy .*= ϵΠx
            fftvy .+= Iy
            fftvy .*= ∂y

            fftvx .+= fftvy
            U[1] .= fftvx

        end

        # Build raw data from physical data.
        # Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
        # Build raw data from physical data.
        # Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
        function mapto(data::InitialData)
            #U = cat(Π⅔.* fft(data.η(x,y)), Π⅔.*fft(data.vx(x,y)), Π⅔.*fft(data.vy(x,y)),dims=3)
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

        # Reconstruct physical variables from raw data
        # Return `(η,vx,vy,x,y)`, where
        # - `η` is the surface deformation;
        # - `vx` is the horizontal velocity field;
        # - `vy` is the horizontal velocity field;
        # - `(x,y)` is the matrix of collocation points
        function mapfro(U)
            real(ifft(U[1])), real(ifft(U[2])), real(ifft(U[3])), x, y
        end


        new(label, f!, mapto, mapfro, info)
    end
end
