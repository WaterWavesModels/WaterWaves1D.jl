export relaxedGreenNaghdi

@doc raw"""
    relaxedGreenNaghdi(param; kwargs...)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
a relaxed Green-Naghdi model proposed by [N. Favrie and S. Gavrilyuk](@cite FavrieGavrilyuk2017) 
or [C. Escalante, M. Dumbser and M. Castro](@cite EscalanteDumbserCastro2019) 
and [G. Richard](@cite Richard2021).

In the latter works the model reads
```math
  \left\{\begin{array}{l}
  ∂_tη+∂_x (hu)=0,\\[1ex]
  h(∂_tu+ϵu∂_x u+∂_x η)+aμ ∂_x (hp)=0,\\[1ex]
  h(∂_tp +ϵu∂_x p)+a(2w+h∂_xu)=0,\\[1ex]
	h(∂_tw+uϵ∂_x w) = a\tfrac{3}{2}p.
  \end{array}\right.
```
where ``h=1 + ϵ η`` is the water depth, ``η`` the surface deformation, ``u`` the horizontal velocity and ``p`` and ``w`` are expected to approximate the layer-veraged pressure and vertical velocity when ``a\gg1``.
Notice the system has four evolution equations, in particular initial data for the augmented variables ``p`` and ``w`` must be suitably chosen. 


# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- the relaxation parameter `a`;
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `FG`: if `true` (default is `false`), compute the Favrie-Gavrilyuk model, otherwise compute the Escalante-Dumbser-Castro model;
- `id`: `∈{0,1,2}` and represent the level of preparation of the initial data (default is `1`);
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `iterative`: solve the elliptic problem (to construct initial data) through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orszag rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Favrie-Gavrilyuk"` if `FG==true`, `"Escalante-Dumbser-Castro"` otherwise);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `relaxedGreenNaghdi.f!` to be called in explicit time-integration solvers;
2. a function `relaxedGreenNaghdi.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `relaxedGreenNaghdi.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
  - `η` is the values of surface deformation at collocation points `x`;
  - `v` is the derivative of the trace of the velocity potential at `x`;
4. additionally, a handy function `relaxedGreenNaghdi.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u,p,w)`, where
  - `u` corresponds to the layer-averaged horizontal velocity.
  - `p` corresponds to the relaxed (artificial) layer-averaged non-hydrostatic pressure;
  - `w` corresponds to the relaxed (artificial) layer-averaged vertical velocity.

"""
mutable struct relaxedGreenNaghdi <: AbstractModel

    label::String
    f!::Function
    mapto::Function
    mapfro::Function
    mapfrofull::Function
    info::String

    function relaxedGreenNaghdi(
            param::NamedTuple; FG = false,
            id = 1,
            mesh = Mesh(param),
            dealias = 0,
            ktol = 0,
            iterate = true,
            gtol = 1.0e-14,
            precond = true,
            restart = nothing,
            maxiter = nothing,
            label = nothing
        )
        # Set up
        a = param.a
        μ = param.μ
        ϵ = param.ϵ

        if isnothing(maxiter)
            maxiter = mesh.N
        end
        if isnothing(restart)
            restart = min(20, mesh.N)
        end
        if isnothing(label)
            if FG == true
                label = "Favrie-Gavrilyuk"
            else
                label = "Escalante-Dumbser-Castro"
            end
        end


        # Print information
        info = "$label model.\n"
        info *= "├─Relaxation parameter a=$a.\n"
        info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ.\n"

        if dealias == 0
            info *= "├─No dealiasing. "
        else
            info *= "├─Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
        end
        if ktol == 0
            info *= "No Krasny filter. "
        else
            info *= "Krasny filter with tolerance $ktol."
        end
        if iterate == true
            if precond == false
                out = "out"
            else
                out = ""
            end
            info *= "\n└─Elliptic problem solved with GMRES method with$out preconditioning, \
            tolerance $gtol, maximal number of iterations $maxiter, restart after $restart iterations \
            (consider `iterate=false` for non-iterative method). "
        else
            info *= "\n└─Elliptic problem solved with standard LU factorization \
            (consider `iterate=true` for faster results). "
        end
        info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

        # Pre-allocate useful data
        k = mesh.k
        x = mesh.x
        x₀ = mesh.x[1]

        ∂ₓ = 1im * k

        if precond == true
            Precond = Diagonal(1 .+ μ / 3 * k .^ 2)
        elseif precond == false
            Precond = Diagonal(ones(size(k)))
        else
            Precond = precond
        end
        if dealias == 0
            Π⅔ = ones(size(k)) # no dealiasing (Π⅔=Id)
        else
            K = (mesh.kmax - mesh.kmin) / (2 + dealias)
            Π⅔ = abs.(k) .<= K # Dealiasing low-pass filter
        end
        FFT = exp.(-1im * k * (x .- x₀)')
        IFFT = exp.(1im * k * (x .- x₀)') / length(x)
        Dx = IFFT * Diagonal(Π⅔ .* ∂ₓ) * FFT
        Id = Diagonal(ones(size(x)))
        h = zeros(Complex{Float64}, mesh.N)
        u, p, w, fftη, fftu, fftv, fftp, fftw, hFG = (similar(h),) .* ones(9)
        L = similar(FFT)


        # Evolution equations are ∂t U = f(U)
        function f!(U; a = a)
            fftη .= U[1]
            h .= 1 .+ ϵ * ifft(fftη)
            fftu .= U[2]; u .= ifft(fftu)
            fftp .= U[3]; p .= ifft(fftp)
            fftw .= U[4]; w .= ifft(fftw)
            if FG == true
                hFG .= h .* (1 .- ϵ * p / a)
            else
                hFG .= h
            end

            U[1] .= -∂ₓ .* Π⅔ .* fft(h .* u)
            U[2] .= -∂ₓ .* Π⅔ .* (fftη .+ ϵ / 2 * fft(u .^ 2)) - a * μ * Π⅔ .* fft(1 ./ h .* ifft(∂ₓ .* fft(hFG .* p)))
            U[3] .= -a * Π⅔ .* fft((2 * w .+ hFG .* ifft(∂ₓ .* fftu)) ./ h) - ϵ * Π⅔ .* fft(u .* ifft(∂ₓ .* fftp))
            U[4] .= 3 / 2 * a * Π⅔ .* fft(p ./ h) - ϵ * Π⅔ .* fft(u .* ifft(∂ₓ .* fftw))
            for u in U
                u[abs.(u) .< ktol] .= 0
            end
            return
        end

        # Build raw data from physical data.
        # Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
        function mapto(data::InitialData)
            fftη .= fft(data.η(x))
            fftv .= fft(data.v(x))
            h .= 1 .+ ϵ * ifft(fftη)
            if iterate == false
                L .= Id - μ / 3 * Diagonal(Π⅔) * FFT * Diagonal(1 ./ h) * Dx * Diagonal(h .^ 3) * IFFT * Diagonal(Π⅔ .* ∂ₓ)
                fftu .= L \ fftv
            else #if iterate == true
                function LL(hatu)
                    return hatu - μ / 3 * Π⅔ .* fft(1 ./ h .* ifft(Π⅔ .* ∂ₓ .* fft(h .^ 3 .* ifft(Π⅔ .* ∂ₓ .* hatu))))
                end
                fftu .= gmres(
                    LinearMap(LL, length(h); issymmetric = false, ismutating = false), fftv;
                    restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol
                )
            end
            if id >= 2
                if iterate == false
                    L .= FFT * Diagonal(3 ./ h .^ 3) * IFFT - μ * Diagonal(Π⅔ .* ∂ₓ) * FFT * Diagonal(1 ./ h) * IFFT * Diagonal(Π⅔ .* ∂ₓ)
                    fftp .= fft(ifft(Π⅔ .* (L \ (Π⅔ .* (∂ₓ .^ 2 .* fftη + 2 * fft(ifft(Π⅔ .* ∂ₓ .* fftu) .^ 2)))) / a) ./ h)
                elseif iterate == true
                    function ll(hatu)
                        return Π⅔ .* fft(3 ./ h .^ 3 .* ifft(hatu)) - μ * Π⅔ .* ∂ₓ .* fft(1 ./ h .* ifft(Π⅔ .* ∂ₓ .* hatu))
                    end
                    fftp .= fft(
                        ifft(
                            Π⅔ .* gmres(
                                LinearMap(ll, length(h); issymmetric = false, ismutating = false), Π⅔ .* (∂ₓ .^ 2 .* fftη + 2 * fft(ifft(Π⅔ .* ∂ₓ .* fftu) .^ 2));
                                restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol
                            ) / a
                        ) ./ h
                    )
                elseif iterate == 1 / 2
                    U₀ = [Π⅔ .* fftη Π⅔ .* fftu fft(zero(x)) Π⅔ .* (-1 / 2 * fft(h .* ifft(∂ₓ .* fftu)))]
                    dt = param.dt .^ 2
                    # Calcul U(dt)
                    U0 = copy(U₀)
                    f!(U0; a = a^2)
                    U1 = copy(U0)

                    U0 .= U₀ .+ dt / 2 .* U1
                    f!(U0; a = a^2)
                    U2 = copy(U0)

                    U0 .= U₀ .+ dt / 2 .* U2
                    f!(U0; a = a^2)
                    U3 = copy(U0)

                    U0 .= U₀ .+ dt .* U3
                    f!(U0; a = a^2)
                    U4 = copy(U0)

                    U₊ = U₀ + dt / 6 .* (U1 + 2 * U2 + 2 * U3 + U4)

                    # Calcul U(-dt)
                    U0 = copy(U₀)
                    f!(U0; a = a^2)
                    U1 = -copy(U0)

                    U0 .= U₀ .+ dt / 2 .* U1
                    f!(U0; a = a^2)
                    U2 = -copy(U0)

                    U0 .= U₀ .+ dt / 2 .* U2
                    f!(U0; a = a^2)
                    U3 = -copy(U0)

                    U0 .= U₀ .+ dt .* U3
                    f!(U0; a = a^2)
                    U4 = -copy(U0)

                    U₋ = U₀ + dt / 6 .* (U1 + 2 * U2 + 2 * U3 + U4)

                    d2th = ifft(U₊[1] + U₋[1] - 2 * U₀[1]) / dt^2
                    dth = ifft(U₊[1] - U₋[1]) / (2 * dt)
                    dtdxh = ifft(∂ₓ .* U₊[1] - ∂ₓ .* U₋[1]) / (2 * dt)
                    d2xh = ifft(∂ₓ .* ∂ₓ .* U₀[1])
                    dtu = ifft(U₊[2] - U₋[2]) / (2 * dt)
                    dxu = ifft(∂ₓ .* U₀[2])
                    dxh = ifft(∂ₓ .* U₀[1])
                    u = ifft(U₀[2])

                    # Need to relaod everything since this has been modified by f!
                    fftη .= fft(data.η(x))
                    fftv .= fft(data.v(x))
                    h .= 1 .+ ϵ * ifft(fftη)
                    fftu .= fft(u)

                    fftp = 1 / 3 / a * Π⅔ .* fft(h .* (d2th + dtu .* dxh + 2 * u .* dtdxh + u .* dxu .* dxh + u .* u .* d2xh))


                end

            else
                fftp .= fft(zero(x))
            end
            if id == 0
                fftw = fft(zero(x))
            else
                if FG == true
                    hFG .= h .* (1 .- ϵ * ifft(fftp))
                else
                    hFG .= h
                end
                fftw .= -1 / 2 * fft(hFG .* ifft(∂ₓ .* fftu))
            end
            U = [Π⅔ .* fftη, Π⅔ .* fftu, Π⅔ .* fftp, Π⅔ .* fftw]
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
            fftη .= U[1]
            h .= 1 .+ ϵ * ifft(fftη)
            fftu .= U[2]
            fftv .= fftu - μ / 3 * Π⅔ .* fft(1 ./ h .* ifft(Π⅔ .* ∂ₓ .* fft(h .^ 3 .* ifft(Π⅔ .* ∂ₓ .* fftu))))
            return real(ifft(fftη)), real(ifft(fftv)), mesh.x
        end
        # Return `(η,v,u,p,w,x)`, where
        # - `η` is the surface deformation;
        # - `v` is the derivative of the trace of the velocity potential;
        # - `u` corresponds to the layer-averaged horizontal velocity.
        # - `p` corresponds to the relaxed (artificial) layer-averaged non-hydrostatic pressure;
        # - `w` corresponds to the relaxed (artificial) layer-averaged vertical velocity.
        # - `x` is the vector of collocation points
        function mapfrofull(U)
            fftη .= U[1]
            h .= 1 .+ ϵ * ifft(fftη)
            fftu .= U[2]
            fftv .= fftu - μ / 3 * Π⅔ .* fft(1 ./ h .* ifft(Π⅔ .* ∂ₓ .* fft(h .^ 3 .* ifft(Π⅔ .* ∂ₓ .* fftu))))
            return real(ifft(fftη)), real(ifft(fftv)), real(ifft(fftu)), real(ifft(U[3])), real(ifft(U[4])), mesh.x
        end

        return new(label, f!, mapto, mapfro, mapfrofull, info)
    end


end
