# loading the packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Pkg.add("Distributions")
# Pkg.add("Plots")
# Pkg.add("LaTeXStrings")
# Pkg.add("Random")
# Pkg.add("Printf")

using Random
using Distributions
using Plots
using LaTeXStrings
using Statistics
using Printf

function plot_kwargs()
    fontsizes = (
        xtickfontsize = 14,
        ytickfontsize = 14,
        xguidefontsize = 16,
        yguidefontsize = 16,
        legendfontsize = 14,
    )
    (;
        linewidth = 3,
        markersize = 6,
        markerstrokewidth = 0,
        fontsizes...,
        size = (600, 500),
    )
end

figdir = joinpath(@__DIR__, "Plots")
mkpath(figdir)

# ─── Parameter ────────────────────────────────────────────────────────────────
N = 100
r0 = 0.05   # minimal distance
rc = 2.5    # Cutoff-Radius LJ
ε = 1.0
σ = 1.0

# a)
function init_particles(N, L, r0)
    data = Tuple{Float64,Float64}[]
    while length(data) < N
        x = rand(Uniform(0, L))
        y = rand(Uniform(0, L))
        if isempty(data)
            push!(data, (x, y))
        else
            valid = true
            for p in data
                if sqrt((x - p[1])^2 + (y - p[2])^2) < r0
                    valid = false
                    break
                end
            end
            valid && push!(data, (x, y))
        end
    end
    return [p[1] for p in data], [p[2] for p in data]
end

function apply_pbc(x, L)
    x < 0   && return x + L
    x >= L  && return x - L
    return x
end

function minimum_image(xij, L)
    xij >  L/2 && return xij - L
    xij < -L/2 && return xij + L
    return xij
end

# b)

function U_LJ(r; ε=1.0, σ=1.0)
    return 4ε * ((σ/r)^12 - (σ/r)^6)
end

function U_shift(r, rc; ε=1.0, σ=1.0)
    r < rc && return U_LJ(r; ε, σ) - U_LJ(rc; ε, σ)
    return 0.0
end

function U_tot(x, y, L; rc=2.5, ε=1.0, σ=1.0)
    U = 0.0
    for i in 1:length(x)
        for j in i+1:length(x)
            xij = minimum_image(x[j] - x[i], L)
            yij = minimum_image(y[j] - y[i], L)
            rij = sqrt(xij^2 + yij^2)
            U += U_shift(rij, rc; ε, σ)
        end
    end
    return U
end

# c)
function energy_change(i, dx, dy, x, y, L; rc=2.5, ε=1.0, σ=1.0)
    E_old = 0.0
    E_new = 0.0
    for j in 1:length(x)
        j == i && continue
        xij_old = minimum_image(x[j] - x[i],      L)
        yij_old = minimum_image(y[j] - y[i],      L)
        xij_new = minimum_image(x[j] - (x[i]+dx), L)
        yij_new = minimum_image(y[j] - (y[i]+dy), L)
        E_old += U_shift(sqrt(xij_old^2 + yij_old^2), rc; ε, σ)
        E_new += U_shift(sqrt(xij_new^2 + yij_new^2), rc; ε, σ)
    end
    return E_new - E_old
end

# d)
function sweep!(x, y, T, Δ, L; rc=2.5, ε=1.0, σ=1.0)
    accepted = 0
    for i in 1:length(x)
        dx = rand(Uniform(-Δ, Δ))
        dy = rand(Uniform(-Δ, Δ))
        dE = energy_change(i, dx, dy, x, y, L; rc, ε, σ)
        if dE < 0 || rand() < exp(-dE / T)
            x[i] = apply_pbc(x[i] + dx, L)
            y[i] = apply_pbc(y[i] + dy, L)
            accepted += 1
        end
    end
    return accepted / length(x)
end

# e)
L  = 10.0
T  = 1.0
Δ  = 0.1
n_eq = 1000
n_samp = 10_000
m  = 10

x_komp, y_komp = init_particles(N, L, r0)

for _ in 1:n_eq
    sweep!(x_komp, y_komp, T, Δ, L)
end

energies = Float64[]
acceptance_ratios = Float64[]

for s in 1:n_samp
    ratio = sweep!(x_komp, y_komp, T, Δ, L)
    push!(acceptance_ratios, ratio)
    mod(s, m) == 0 && push!(energies, U_tot(x_komp, y_komp, L))
end

println("⟨E⟩/N = ", mean(energies) / N)
println("Standardfehler = ", std(energies) / sqrt(length(energies)))
println("Akzeptanzrate  = ", mean(acceptance_ratios))

# Delta-Vergleich
println("\n── Delta-Vergleich ──")
for Δ_test in [0.001, 0.01, 0.1, 0.5]
    x_t, y_t = copy(x_komp), copy(y_komp)
    for _ in 1:n_eq
        sweep!(x_t, y_t, T, Δ_test, L)
    end
    E_tmp  = Float64[]
    ar_tmp = Float64[]
    for s in 1:n_samp
        r = sweep!(x_t, y_t, T, Δ_test, L)
        push!(ar_tmp, r)
        mod(s, m) == 0 && push!(E_tmp, U_tot(x_t, y_t, L))
    end
    @printf("Δ = %.3f | ⟨E⟩/N = %.4f | SE = %.4f | acc = %.4f\n",
        Δ_test, mean(E_tmp)/N, std(E_tmp)/sqrt(length(E_tmp)), mean(ar_tmp))
end

# f)
function virial(x, y, L; rc=2.5, ε=1.0, σ=1.0)
    W = 0.0
    for i in 1:length(x)
        for j in i+1:length(x)
            xij = minimum_image(x[j] - x[i], L)
            yij = minimum_image(y[j] - y[i], L)
            rij = sqrt(xij^2 + yij^2)
            if rij < rc
                # f·r = dU/dr * r = 4ε(12(σ/r)^12 - 6(σ/r)^6)
                W += 4ε * (12*(σ/rij)^12 - 6*(σ/rij)^6)
            end
        end
    end
    return 0.5 * W
end

function pressure(rho, T, x, y, L; rc=2.5, ε=1.0, σ=1.0)
    A = L^2
    W = virial(x, y, L; rc, ε, σ)
    ΔP_tail = π * rho^2 * ε * (24/11 * (σ/rc)^12 - 12/5 * (σ/rc)^6)
    return rho * T + W / A + ΔP_tail
end

# Plot P vs rho
T_plot = 2.0
Δ_plot = 0.65
n_eq_p = 1000
n_samp_p = 1000
rho_vals = 0.01:0.02:0.5

pressures_mc = Float64[]

for rho in rho_vals
    L_rho = sqrt(N / rho)
    x_r, y_r = init_particles(N, L_rho, r0)

    for _ in 1:n_eq_p
        sweep!(x_r, y_r, T_plot, Δ_plot, L_rho)
    end

    P_samples = Float64[]
    for s in 1:n_samp_p
        sweep!(x_r, y_r, T_plot, Δ_plot, L_rho)
        if mod(s, m) == 0
            push!(P_samples, pressure(rho, T_plot, x_r, y_r, L_rho))
        end
    end
    push!(pressures_mc, mean(P_samples))
    println("rho = $rho | P = $(pressures_mc[end])")
end

rho_plot = collect(rho_vals)
plt = plot(rho_plot, pressures_mc;
    label="LJ MC", xlabel=L"\rho", ylabel="P", plot_kwargs()...)
plot!(plt, rho_plot, rho_plot .* T_plot;
    label="Ideal gas", linestyle=:dash, plot_kwargs()...)
savefig(plt, joinpath(figdir, "P_vs_rho.pdf"))
