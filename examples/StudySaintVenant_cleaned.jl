# =============================================================================
# StudySaintVenant.jl
# Reproduces the figures in the work of V. Duchêne and J. Marstrander
# on the numerical discretization of quasilinear systems.
# This code has undergone postprocessing to improve formatting, indentation,
# and in-line documentation with the assistance of OpenAI's ChatGPT.
# The scientific content and numerical implementation remain entirely
# the work of the original author(s), V. Duchêne and J. Marstrander, 
# and have not been altered.
# =============================================================================

export IntegrateSV, IntegrateSV2D, Convergence, Convergence2D,
       Figure1, Figure2, Figure3, Figure4, Figure5, Figure6, Figure7

using WaterWaves1D, FFTW, Plots, LinearAlgebra, ProgressMeter

# Set default plot font
plot_font = "Computer Modern"
default(fontfamily = plot_font)

# ------------------- 
# --- 1D PROBLEMS ---
# -------------------

# ----------------------------------------------------------------
# IntegrateSV: Numerical integration of the 1D Saint-Venant system
# ----------------------------------------------------------------

"""
    IntegrateSV(; init, args...)

Integrate the Saint-Venant (shallow water) system in 1D with initial data depending on `init`.

- `init = 1`: η₀(x) = ½·exp(-|x|^α)·exp(-4x²), v₀(x) = 0
- `init = 2`: η₀(x) = (h₀-1)·cos(x), v₀(x) = v₀·sin(x) + sin(N·x)/N²
- `init = 3`: η₀(x) = ½·exp(-4x²), v₀(x) = 2·√(1+ε·η₀) - 2

Keyword arguments:
- `α`: regularity exponent (default 3/2)
- `N`, `h₀`, `v₀`: frequency, minimum depth, velocity for init=2
- `ϵ`: nonlinearity parameter (default 1)
- `L`: half-length of mesh (default π)
- `M`: number of collocation points (default 2⁹)
- `T`: final integration time (default 0.5)
- `dt`: time step (default 1e-5)
- `dealias`: whether to use dealiasing (default true)
- `smooth`: smooth low-pass filter (default false)
- `Ns`: number of stored output times
- `label`: descriptive string

Returns the `problem` object containing the solution.
"""
function IntegrateSV(; init, α=1.5, N=nothing, h₀=1/2, v₀=2,
                      ϵ=1, L=π, M=2^9, T=0.5, dt=1e-5,
                      dealias=true, smooth=false, Ns=nothing,
                      label="Saint-Venant")

    # Set simulation parameters
    param = isnothing(Ns) ?
        (ϵ=ϵ, N=M, L=L, T=T, dt=dt) :
        (ϵ=ϵ, N=M, L=L, T=T, dt=dt, Ns=Ns)

    mesh = Mesh(param)

    # Select initial condition based on `init`
    if init == 1
        init = Init(x -> 0.5 * exp.(-abs.(x).^α) .* exp.(-4x.^2),
                    x -> 0 * cos.(x))
    elseif init == 2
        if isnothing(N)
            N = M ÷ 3
        end
        init = Init(x -> (h₀ - 1) * cos.(x),
                    x -> v₀ .* sin.(x) .+ sin.(N * x) / N^2)
    elseif init == 3
        f(x) = 0.5 * exp.(-4 * x.^2)
        init = Init(x -> f(x),
                    x -> 2 * sqrt.(1 .+ ϵ * f(x)) .- 2)
    else
        @error "argument init must be 1, 2, or 3"
    end

    # Define model, solver and set initial-value problem
    model  = SaintVenant_fast(param; dealias=dealias, smooth=smooth, label=label)
    solver = RK4(model.mapto(init))
    problem = Problem(model, init, param; solver=solver)

    # Run simulation
    solve!(problem)

    # --- Check energy preservation ---

    e(η, v) = 0.5 * (η.^2 .+ (1 .+ ϵ * η) .* v.^2)
    e(η1, v1, η2, v2) = 0.5 * ((η1 - η2) .* (η1 + η2) .+
                               (1 .+ ϵ * η2) .* (v1 - v2) .* (v1 + v2) .+
                               ϵ * (η1 - η2) .* v1.^2)

    ηfin, vfin = solution(problem)
    η0, v0 = init.η(mesh.x), init.v(mesh.x)

    error_energy = abs(sum(e(ηfin, vfin, η0, v0)) / sum(e(η0, v0)))
    @info "normalized preservation of the total energy: $error_energy\\n"

    return problem
end

# ------------------------------------------------------------
# Convergence: Convergence test for the 1D Saint-Venant system
# ------------------------------------------------------------

"""
    Convergence(init; smooth=false, T=0.5)

Compute the convergence rate by solving the Saint-Venant system for
various resolutions and comparing to a reference solution.

Arguments:
- `init = 1`: heap of water (Figures 1 and 2)
- `init = 2`: high-frequency mode (Figure 3)
- `smooth`: use smooth low-pass filter (default: false)
- `T`: final integration time (default: 0.5)

Returns:
- reference problem
- L² relative error array
- H¹ relative error array
"""
function Convergence(; init, smooth=false, name=nothing, T=0.5)
    E0 = Float64[]  # L² errors,
    E1 = Float64[]  # H¹ errors

    if init == 1
        # Reference solution
        reference_problem = IntegrateSV(init=1, α=1.5, ϵ=1, L=π, M=2^15,
                                        T=T, dt=1e-5, dealias=true,
                                        smooth=false, Ns=100, label="reference")
        # Fourier coefficients
        Uref = reference_problem.data.U[end] / 2^15

        for n in 14:-1:6
            # Lower-resolution run
            p = IntegrateSV(init=1, α=1.5, ϵ=1, L=π, M=2^n, T=T, dt=1e-5,
                            dealias=true, smooth=smooth, Ns=1, label="2M=2^$n")
            # Fourier coefficients
            U = [p.data.U[end][1]; p.data.U[end][2]] / 2^n

            Nmodes = 2^(n-1)
            Ucomp = [Uref[1][1:Nmodes]; Uref[1][end-Nmodes+1:end];
                     Uref[2][1:Nmodes]; Uref[2][end-Nmodes+1:end]]

            # L² error
            push!(E0, norm(U - Ucomp, 2) / norm(Ucomp, 2))

            # H¹ error
            k = [0:Nmodes-1; -Nmodes:-1]
            K = [sqrt.(1 .+ k.^2); sqrt.(1 .+ k.^2)]
            push!(E1, norm(K .* U - K .* Ucomp, 2) / norm(K .* Ucomp, 2))
        end
    elseif init == 2
        for n in 14:-1:6
            # Reference solution
            reference_problem = IntegrateSV(init=2, N=(2^n ÷ 3), h₀=0.5, v₀=2,
                                            ϵ=1, L=π, M=2^15, T=T, dt=1e-5,
                                            dealias=true, smooth=smooth,
                                            Ns=1, label="reference")
            # Lower-resolution run
            p = IntegrateSV(init=2, h₀=0.5, v₀=2, ϵ=1, L=π, M=2^n, T=T,
                            dt=1e-5, dealias=true, smooth=smooth,
                            Ns=1, label="2M=2^$n")
            # Fourier coefficients
            Uref = reference_problem.data.U[end] / 2^15
            U = [p.data.U[end][1]; p.data.U[end][2]] / 2^n

            Nmodes = 2^(n-1)
            Ucomp = [Uref[1][1:Nmodes]; Uref[1][end-Nmodes+1:end];
                     Uref[2][1:Nmodes]; Uref[2][end-Nmodes+1:end]]

            # L² error
            push!(E0, norm(U - Ucomp, 2) / norm(Ucomp, 2))
            
            # H¹ error
            k = [0:Nmodes-1; -Nmodes:-1]
            K = [sqrt.(1 .+ k.^2); sqrt.(1 .+ k.^2)]
            push!(E1, norm(K .* U - K .* Ucomp, 2) / norm(K .* Ucomp, 2))
        end
    end

    return reference_problem, E0, E1
end

# ----------------------------------------------------------------------------
# Figure 1: Initial data : heap of water, physical and spectral representation
# ----------------------------------------------------------------------------
function Figure1()
    # Build initial data
    reference_problem, E0, E1 = Convergence(init=1; smooth=false, T=0)
    η0, v0, x = solution(reference_problem,T=0)

    Fig1a = plot(x, [η0 v0], label=["\$\\eta^0\$" "\$u^0\$"], xlabel="\$x\$")
    savefig(Fig1a, "Fig1a.pdf"); savefig(Fig1a, "Fig1a.svg")

    Fig1b = plot(; xlabel="\$2M\$", axis=:log)
    n = 14:-1:6
    M = 2 .^ n
    scatter!(Fig1b, M, E1, label="\$E_1\$", color=1)
    scatter!(Fig1b, M, E0, label="\$E_0\$", color=2)
    plot!(Fig1b, M, M.^(-1), label="", color=1)
    plot!(Fig1b, M, M.^(-2), label="", color=2)
    savefig(Fig1b, "Fig1b.pdf"); savefig(Fig1b, "Fig1b.svg")
end

# --------------------------------------------------------
# Figure 2: Convergence with sharp vs smooth low-pass filter
# --------------------------------------------------------
function Figure2()
    # Sharp low-pass filter
    _, E0, E1 = Convergence(init=1; smooth=false)

    # Smooth low-pass filter
    _, E0_smooth, E1_smooth = Convergence(init=1; smooth=true)

    # Convergence plot
    Fig2 = plot(; xlabel="\$2M\$", axis=:log)
    n = 14:-1:6
    M = 2 .^ n

    scatter!(Fig2, M, E1, label="\$E_1\$, sharp low-pass filter", color=1)
    scatter!(Fig2, M, E0, label="\$E_0\$, sharp low-pass filter", color=2)
    scatter!(Fig2, M, E1_smooth, label="\$E_1\$, smooth low-pass filter", color=3)
    scatter!(Fig2, M, E0_smooth, label="\$E_0\$, smooth low-pass filter", color=4)
    plot!(Fig2, M, M.^(-1), label="", color=1)
    plot!(Fig2, M, M.^(-2), label="", color=2)

    savefig(Fig2, "Fig2.pdf"); savefig(Fig2, "Fig2.svg")

    return round.(log2.(E0[1:end-1] ./ E0[2:end]), digits=2),
           round.(log2.(E1[1:end-1] ./ E1[2:end]), digits=2),
           round.(log2.(E0_smooth[1:end-1] ./ E0_smooth[2:end]), digits=2),
           round.(log2.(E1_smooth[1:end-1] ./ E1_smooth[2:end]), digits=2)
end

# ----------------------------------------------------------
# Figure 3: High-frequency mode initial data and convergence
# ----------------------------------------------------------
function Figure3()
    # Build and plot initial data with high-frequency component
    reference_problem, _, _ = Convergence(init=2; smooth=false, T=0)
    η0, v0, x = reference_problem.model.mapfro(reference_problem.data.U[1])

    Fig3a = plot(x, [η0 v0], label=["\$\\eta^0\$" "\$u^0\$"], xlabel="\$x\$")
    plot!(x, v0 .- 2*sin.(x),
          xlim=(0, π),
          frame=:box,
          color=2,
          inset=bbox(0.6, 0.6, 0.35, 0.25),
          subplot=2,
          label="\$u^0 - 2\\sin(x)\$")
    savefig(Fig3a, "Fig3a.pdf")
    savefig(Fig3a, "Fig3a.svg")

    # Sharp low-pass filter
    _, E0, E1 = Convergence(init=2; smooth=false)

    # Smooth low-pass filter
    _, E0_smooth, E1_smooth = Convergence(init=2; smooth=true)

    # Convergence plot
    Fig3b = plot(; xlabel="\$2M\$", axis=:log, legend=:bottomleft)
    M = 2 .^ (14:-1:6)
    scatter!(Fig3b, M, E1, label="\$E_1\$, sharp", color=1)
    scatter!(Fig3b, M, E0, label="\$E_0\$, sharp", color=2)
    scatter!(Fig3b, M, E1_smooth, label="\$E_1\$, smooth", color=3)
    scatter!(Fig3b, M, E0_smooth, label="\$E_0\$, smooth", color=4)
    plot!(Fig3b, M, M.^(-1), label="", color=1)
    plot!(Fig3b, M, M.^(-2), label="", color=2)

    savefig(Fig3b, "Fig3b.pdf")
    savefig(Fig3b, "Fig3b.svg")
end

# ----------------------------------------------------------
# Figure 4: Instability in simulations with vanishing depth
# ----------------------------------------------------------
function Figure4()
    # Simulations with M = 2^10
    # Smooth low-pass filter
    pb_smooth_M10 = IntegrateSV(init=2, h₀=0, v₀=2, ϵ=1, L=π, M=2^10, T=0.1,
                                dt=1e-5, dealias=true, smooth=true, Ns=1, label="smooth")
    # Sharp low-pass filter
    pb_sharp_M10 = IntegrateSV(init=2, h₀=0, v₀=2, ϵ=1, L=π, M=2^10, T=0.1,
                               dt=1e-5, dealias=true, smooth=false, Ns=1, label="sharp")

    
    # Plot second derivative of the velocity
    η_M10, v_M10, x = solution(pb_sharp_M10)
    ηs_M10, vs_M10, _ = solution(pb_smooth_M10)

    ∂ = 1im * Mesh(x).k
    ∂2v_M10 = real.(ifft(∂.^2 .* fft(v_M10)))
    ∂2vs_M10 = real.(ifft(∂.^2 .* fft(vs_M10)))

    Fig4a = plot(x, ∂2v_M10, label="sharp")
    plot!(x, ∂2vs_M10, label="smooth")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{10}\$")
    xlabel!("\$x\$")
    savefig(Fig4a, "Fig4a.pdf")
    savefig(Fig4a, "Fig4a.svg")

    # Simulations with M = 2^12
    # Smooth low-pass filter
    pb_smooth_M12 = IntegrateSV(init=2, h₀=0, v₀=2, ϵ=1, L=π, M=2^12, T=0.1,
                                dt=1e-5, dealias=true, smooth=true, Ns=1, label="smooth")
    # Sharp low-pass filter
    pb_sharp_M12 = IntegrateSV(init=2, h₀=0, v₀=2, ϵ=1, L=π, M=2^12, T=0.1,
                               dt=1e-5, dealias=true, smooth=false, Ns=1, label="sharp")

    # Plot second derivative of the velocity
    η_M12, v_M12, x = solution(pb_sharp_M12)
    ηs_M12, vs_M12, _ = solution(pb_smooth_M12)

    ∂ = 1im * Mesh(x).k
    ∂2v_M12 = real.(ifft(∂.^2 .* fft(v_M12)))
    ∂2vs_M12 = real.(ifft(∂.^2 .* fft(vs_M12)))

    Fig4b = plot(x, ∂2v_M12, label="sharp")
    plot!(x, ∂2vs_M12, label="smooth")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{12}\$")
    xlabel!("\$x\$")
    savefig(Fig4b, "Fig4b.pdf")
    savefig(Fig4b, "Fig4b.svg")
end

# -----------------------------------------------------------------
# FigureX: Simulations of smooth simple waves after wavebreaking 
# (not appearingin the manuscript)
# -----------------------------------------------------------------
function FigureX()
    # Simulations with M = 2^10 modes
    pb_sharp_M10 = IntegrateSV(init=3, ϵ=1, L=π, M=2^10, T=1, dt=1e-5,
                               dealias=true, smooth=false, Ns=100, label="sharp")
    pb_smooth_M10 = IntegrateSV(init=3, ϵ=1, L=π, M=2^10, T=1, dt=1e-5,
                                dealias=true, smooth=true, Ns=100, label="smooth")

    plot([pb_sharp_M10, pb_smooth_M10], T=1)

    # Simulations with M = 2^12 modes
    pb_sharp_M12 = IntegrateSV(init=3, ϵ=1, L=π, M=2^12, T=1, dt=1e-5,
                               dealias=true, smooth=false, Ns=100, label="sharp")
    pb_smooth_M12 = IntegrateSV(init=3, ϵ=1, L=π, M=2^12, T=1, dt=1e-5,
                                dealias=true, smooth=true, Ns=100, label="smooth")

    plot([pb_sharp_M12, pb_smooth_M12], T=1)
end

# -------------------
# --- 2D PROBLEMS ---
# -------------------

# ---------------------------------------------------------------
# IntegrateSV2D: Numerical integration of the 2D Saint-Venant system
# ---------------------------------------------------------------
"""
    IntegrateSV2D(; args...)

Integrate in time the 2D Saint-Venant (shallow water) system.

Initial data is
- ζ(x, y) = (h₀-1)·cos(x)·cos.(y)'
- u(x, y) = u₀·sin(x)·cos.(y)+ u₁·sin(N·x)·cos(N·y)/N^s
- v(x, y) = v₀·cos(x)·sin.(y)+ v₁·cos(N·x)·sin.(N·y)/N^s

Keyword arguments:
- `N, s, u₀, v₀, u₁, v₁`: Initial data params (frequencies and amplitudes)
- `ϵ`: Nonlinearity parameter (default 1)
- `L`: Half-length of domain (default π)
- `M`: Number of collocation points (default 2⁹)
- `T`: Final time (default 0.1)
- `dt`: Timestep (default 5e-5)
- `dealias`: Use dealiasing (default true)
- `smooth`: Use smooth spectral filter (default false)
- `hamiltonian`: Use Hamiltonian system (default false)
- `Ns`: Number of output time steps (default 1)
- `label`: String label for identification

Returns:
- `problem`: The `Problem` object containing the simulation result
"""
function IntegrateSV2D(; N=nothing, s=2, h₀=1/2,
                        u₀=1/2, v₀=-1/2, u₁=1, v₁=1,
                        ϵ=1, L=π, M=2^9, T=0.1, dt=5e-5,
                        dealias=true, smooth=false,
                        hamiltonian=false, Ns=1,
                        label="Saint-Venant")

    # Define parameters
    param = isnothing(Ns) ?
        (ϵ=ϵ, N=M, L=L, T=T, dt=dt) :
        (ϵ=ϵ, N=M, L=L, T=T, dt=dt, Ns=Ns)

    mesh = Mesh(param)
    if isnothing(N)
        N = M ÷ 3
    end

    # Define initial conditions
    ζ(x, y)   = (h₀ - 1) * cos.(x) .* cos.(y)'
    ux(x, y)  = u₀ * cos.(y') .* sin.(x) .+ u₁ * cos.(N * y') .* sin.(N * x) / N^s
    uy(x, y)  = v₀ * sin.(y') .* cos.(x) .+ v₁ * sin.(N * y') .* cos.(N * x) / N^s

    init = Init2D(ζ, ux, uy)

    # Define model, solver and set initial-value problem
    model = SaintVenant2D_fast(param; dealias=dealias, hamiltonian=hamiltonian, smooth=smooth)
    solver = RK4(model.mapto(init))
    problem = Problem(model, init, param; solver=solver)

    # Run simulation
    solve!(problem)

    # --- Energy check ---
    e(η, vx, vy) = 0.5 * (η.^2 .+ (1 .+ ϵ * η) .* (vx.^2 .+ vy.^2))
    e(η1, vx1, vy1, η2, vx2, vy2) = 0.5 * (
        (η1 - η2) .* (η1 + η2) .+
        (1 .+ ϵ * η2) .* (vx1 - vx2) .* (vx1 + vx2) .+
        (1 .+ ϵ * η2) .* (vy1 - vy2) .* (vy1 + vy2) .+
        ϵ * (η1 - η2) .* vx1.^2 .+ ϵ * (η1 - η2) .* vy1.^2
    )

    ηfin, vxfin, vyfin = solution(problem)
    η0 = init.η(mesh.x, mesh.x)
    vx0 = init.vx(mesh.x, mesh.x)
    vy0 = init.vy(mesh.x, mesh.x)

    error_energy = abs(sum(e(ηfin, vxfin, vyfin, η0, vx0, vy0)) / sum(e(η0, vx0, vy0)))
    @info "normalized preservation of the total energy: $error_energy\\n"

    return problem
end

# --------------------------------------------------------------
# Convergence2D: Convergence test for the 2D Saint-Venant system
# --------------------------------------------------------------
"""
    Convergence2D(; smooth=false, hamiltonian=false, T=0.1)

Compute convergence rates by comparing lower-resolution solutions
to a high-resolution reference for 2D Saint-Venant system.

Returns:
- L² error array
- H¹ error array
"""
function Convergence2D(; smooth=false, hamiltonian=false, name=nothing, T=0.1)
    E0 = Float64[]  # L² errors
    E1 = Float64[]  # H¹ errors

    for n in 9:-1:6
        # Reference solution
        reference_problem = IntegrateSV2D(M=2^10, N=(2^n)÷3, T=T,
                                          smooth=smooth, hamiltonian=hamiltonian,
                                          Ns=1, label="reference") 
        # Fourier coefficients
        Uref = reference_problem.data.U[end] / 2^20
        Uref1 = fftshift(Uref[1])
        Uref2 = fftshift(Uref[2])
        Uref3 = fftshift(Uref[3])
        
        # Lower-resolution run
        p = IntegrateSV2D(M=2^n, N=(2^n)÷3, T=T,
                          smooth=smooth, hamiltonian=hamiltonian,
                          Ns=1, label="2M=2^$n")
        # Fourier coefficients
        U = p.data.U[end] / 2^(2n)
        U1 = fftshift(U[1])
        U2 = fftshift(U[2])
        U3 = fftshift(U[3])

        Nmodes = 2^(n-1)
        Nref = 2^9
        Ucomp1 = Uref1[Nref-Nmodes+1:Nref+Nmodes, Nref-Nmodes+1:Nref+Nmodes]
        Ucomp2 = Uref2[Nref-Nmodes+1:Nref+Nmodes, Nref-Nmodes+1:Nref+Nmodes]
        Ucomp3 = Uref3[Nref-Nmodes+1:Nref+Nmodes, Nref-Nmodes+1:Nref+Nmodes]

        # L² error
        push!(E0, norm([norm(U1 - Ucomp1); norm(U2 - Ucomp2); norm(U3 - Ucomp3)]) /
                  norm([norm(Ucomp1); norm(Ucomp2); norm(Ucomp3)]))

        # H¹ error
        kx = range(-Nmodes, Nmodes-1)
        ky = kx'
        k = sqrt.(1 .+ kx.^2 .+ ky.^2)
        push!(E1, norm([norm(k .* U1 - k .* Ucomp1);
                        norm(k .* U2 - k .* Ucomp2);
                        norm(k .* U3 - k .* Ucomp3)]) /
                  norm([norm(k .* Ucomp1);
                        norm(k .* Ucomp2);
                        norm(k .* Ucomp3)]))
    end

    return E0, E1
end

# ----------------------------------------------------------------
# Figure 5: 2D convergence (sharp and smooth, Hamiltonian and not)
# ----------------------------------------------------------------
function Figure5()
    # Non-Hamiltonian system, sharp and smooth low-pass filters
    E0, E1 = Convergence2D(T=0.1, smooth=false, hamiltonian=false)
    E0s, E1s = Convergence2D(T=0.1, smooth=true, hamiltonian=false)

    Fig5a = plot(; xlabel="\$M\$", axis=:log, legend=:bottomleft)
    M = 2 .^ (9:-1:6)
    scatter!(Fig5a, M, E1, label="\$E_1\$, sharp low-pass filter", color=1)
    scatter!(Fig5a, M, E0, label="\$E_0\$, sharp low-pass filter", color=2)
    scatter!(Fig5a, M, E1s, label="\$E_1\$, smooth low-pass filter", color=3)
    scatter!(Fig5a, M, E0s, label="\$E_0\$, smooth low-pass filter", color=4)
    plot!(Fig5a, M, 4 .* M.^(-1), label="", color=1)
    plot!(Fig5a, M, 12 .* M.^(-2), label="", color=2)
    savefig(Fig5a, "Fig5a.pdf"); savefig(Fig5a, "Fig5a.svg")

    # Hamiltonian system, sharp and smooth low-pass filters
    E0, E1 = Convergence2D(T=0.1, smooth=false, hamiltonian=true)
    E0s, E1s = Convergence2D(T=0.1, smooth=true, hamiltonian=true)

    Fig5b = plot(; xlabel="\$M\$", axis=:log, legend=:bottomleft)
    scatter!(Fig5b, M, E1, label="\$E_1\$, sharp", color=1)
    scatter!(Fig5b, M, E0, label="\$E_0\$, sharp", color=2)
    scatter!(Fig5b, M, E1s, label="\$E_1\$, smooth", color=3)
    scatter!(Fig5b, M, E0s, label="\$E_0\$, smooth", color=4)
    plot!(Fig5b, M, 4 .* M.^(-1), label="", color=1)
    plot!(Fig5b, M, 12 .* M.^(-2), label="", color=2)
    savefig(Fig5b, "Fig5b.pdf"); savefig(Fig5b, "Fig5b.svg")
end

# ------------------------------------------------------
# Figure 6: Instability for 2D simulation with h₀ < 0
# ------------------------------------------------------
function Figure6()
    # M = 2^9
    # Smooth low-pass filter
    problem = IntegrateSV2D(h₀=-0.1, u₀=0.5, v₀=-0.5, u₁=1, v₁=1,
                            M=2^9, T=0.1, dt=5e-5, smooth=false,
                            hamiltonian=false, Ns=1)
    η, vx, vy, x, y = solution(problem)

    # Sharp low-pass filter
    problem_smooth = IntegrateSV2D(h₀=-0.1, u₀=0.5, v₀=-0.5, u₁=1, v₁=1,
                                   M=2^9, T=0.1, dt=5e-5, smooth=true,
                                   hamiltonian=false, Ns=1)

    ηs, vxs, vys, = solution(problem_smooth)

    # Plot second derivative of the velocity
    mesh = Mesh(x)
    kx, ky = mesh.k, mesh.k
    ∂x, ∂y = 1im * kx, 1im * ky'

    ∂fx(f) = real.(ifft(∂x .* fft(f, 1), 1))
    ∂fy(f) = real.(ifft(∂y .* fft(f, 2), 2))

    Fig6a = plot(x, y, ∂fx(∂fx(vx)), st=:surface, label="sharp")
    plot!(xlabel="x", ylabel="y")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{9}\$")
    savefig(Fig6a, "Fig6a.pdf"); savefig(Fig6a, "Fig6a.svg")

    Fig6b = plot(x, y, ∂fx(∂fx(vxs)), st=:surface, label="smooth")
    plot!(xlabel="x", ylabel="y")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{9}\$")
    savefig(Fig6b, "Fig6b.pdf"); savefig(Fig6b, "Fig6b.svg")

    # M = 2^10...
    # Smooth low-pass filter
    problem = IntegrateSV2D(h₀=-0.1, u₀=0.5, v₀=-0.5, u₁=1, v₁=1,
                            M=2^10, T=0.1, dt=5e-5, smooth=false,
                            hamiltonian=false, Ns=1)
    η, vx, vy, x, y = solution(problem)
    # Sharp low-pass filter
    problem_smooth = IntegrateSV2D(h₀=-0.1, u₀=0.5, v₀=-0.5, u₁=1, v₁=1,
                                   M=2^10, T=0.1, dt=5e-5, smooth=true,
                                   hamiltonian=false, Ns=1)
    ηs, vxs, vys, = solution(problem_smooth)

    # Plot second derivative of the velocity
    mesh = Mesh(x)
    kx, ky = mesh.k, mesh.k
    ∂x, ∂y = 1im * kx, 1im * ky'

    ∂fx(f) = real.(ifft(∂x .* fft(f, 1), 1))
    ∂fy(f) = real.(ifft(∂y .* fft(f, 2), 2))

    Fig6c = plot(x, y, ∂fx(∂fx(vx)), st=:surface, label="sharp")
    plot!(xlabel="x", ylabel="y")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{10}\$")
    savefig(Fig6c, "Fig6c.pdf"); savefig(Fig6c, "Fig6c.svg")

    Fig6d = plot(x, y, ∂fx(∂fx(vxs)), st=:surface, label="smooth")
    plot!(xlabel="x", ylabel="y")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{10}\$")
    savefig(Fig6d, "Fig6d.pdf"); savefig(Fig6d, "Fig6d.svg")end

# -------------------------------------------------------------
# Figure 7: Instability for Hamiltonian system simulation in 2D
# -------------------------------------------------------------
function Figure7()
    # M = 2^9
    # Smooth low-pass filter
    problem = IntegrateSV2D(h₀=0.5, u₀=2, v₀=-2, u₁=1, v₁=-1,
                            M=2^9, T=0.1, dt=5e-5, smooth=false,
                            hamiltonian=true, Ns=1)
    η, vx, vy, x, y = solution(problem)

    # Sharp low-pass filter
    problem_smooth = IntegrateSV2D(h₀=0.5, u₀=2, v₀=-2, u₁=1, v₁=-1,
                                   M=2^9, T=0.1, dt=5e-5, smooth=true,
                                   hamiltonian=true, Ns=1)
    ηs, vxs, vys, = solution(problem_smooth)

    # Plot second derivative of the velocity
    mesh = Mesh(x)
    kx, ky = mesh.k, mesh.k
    ∂x, ∂y = 1im * kx, 1im * ky'

    ∂fx(f) = real.(ifft(∂x .* fft(f, 1), 1))
    ∂fy(f) = real.(ifft(∂y .* fft(f, 2), 2))

    Fig7a = plot(x, y, ∂fx(∂fx(vx)), st=:surface, label="sharp")
    plot!(xlabel="x", ylabel="y")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{9}\$")
    savefig(Fig7a, "Fig7a.pdf"); savefig(Fig7a, "Fig7a.svg")

    Fig7b = plot(x, y, ∂fx(∂fx(vxs)), st=:surface, label="smooth")
    plot!(xlabel="x", ylabel="y")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{9}\$")
    savefig(Fig7b, "Fig7b.pdf"); savefig(Fig7b, "Fig7b.svg")

    # M = 2^10
    # Smooth low-pass filter
    problem = IntegrateSV2D(h₀=0.5, u₀=2, v₀=-2, u₁=1, v₁=-1,
                            M=2^10, T=0.1, dt=5e-5, smooth=false,
                            hamiltonian=true, Ns=1)
    η, vx, vy, x, y = solution(problem)

    # Sharp low-pass filter
    problem_smooth = IntegrateSV2D(h₀=0.5, u₀=2, v₀=-2, u₁=1, v₁=-1,
                                   M=2^10, T=0.1, dt=5e-5, smooth=true,
                                   hamiltonian=true, Ns=1)
    ηs, vxs, vys, = solution(problem_smooth)

    # Plot second derivative of the velocity
    mesh = Mesh(x)
    kx, ky = mesh.k, mesh.k
    ∂x, ∂y = 1im * kx, 1im * ky'

    ∂fx(f) = real.(ifft(∂x .* fft(f, 1), 1))
    ∂fy(f) = real.(ifft(∂y .* fft(f, 2), 2))

    Fig7c = plot(x, y, ∂fx(∂fx(vx)), st=:surface, label="sharp")
    plot!(xlabel="x", ylabel="y")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{10}\$")
    savefig(Fig7c, "Fig7c.pdf"); savefig(Fig7c, "Fig7c.svg")

    Fig7d = plot(x, y, ∂fx(∂fx(vxs)), st=:surface, label="smooth")
    plot!(xlabel="x", ylabel="y")
    title!("\$\\partial_x^2 u_N, \\quad 2M = 2^{10}\$")
    savefig(Fig7d, "Fig7d.pdf"); savefig(Fig7d, "Fig7d.svg")
end

nothing

