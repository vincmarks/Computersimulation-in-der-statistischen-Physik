 # wsl, vmd
 
# loading the packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Pkg.add("Distributions")
# Pkg.add("Plots")
# Pkg.add("LaTeXStrings")
# Pkg.add("Random")
# Pkg.add("DelimitedFiles")
# Pkg.add("Statistics")

using Random
using Distributions
using Plots
using LaTeXStrings
using DelimitedFiles
using Statistics

function plot_kwargs()
        fontsizes = (
            xtickfontsize = 14,
            ytickfontsize = 14,
            xguidefontsize = 16,
            yguidefontsize = 16,
            legendfontsize = 14,
        )
        (;
            #gridlinewidth = 2,            
            linewidth = 3,
            markersize = 6,
            markerstrokewidth = 0,
            fontsizes...,
            size = (600, 500),
        )
    end

figdir = joinpath(@__DIR__, "Plots")

############################################

# Exercise 2: Scaling behavior of polymer chains: R²_ee ~ N^(2ν)

datadir = joinpath(@__DIR__, "CODE_PIVOT", "TEST")
N_vals_mc = [200, 400, 600, 800, 1000]

R2_means = Float64[]
R2_errs  = Float64[]

for N in N_vals_mc
    data = readdlm(joinpath(datadir, "analyse_N$N"))
    col2 = data[:, 2]  # R²_ee is the second column, as explained in exercise
    push!(R2_means, mean(col2))
    push!(R2_errs,  std(col2) / sqrt(length(col2)))
    println("N=$N:  <R2_ee> = $(round(mean(col2), digits=2)) +/- $(round(std(col2)/sqrt(length(col2)), digits=2))")
end

# power-law fit in log-log: log(R²) = log(A) + 2ν·log(N)
log_N  = log.(Float64.(N_vals_mc))
log_R2 = log.(R2_means)
X      = [ones(length(log_N)) log_N]
coeffs = X \ log_R2
intercept, slope = coeffs[1], coeffs[2]
nu_fit = slope / 2

N_fit   = range(150, 1100, length=200)
R2_fit  = exp(intercept) .* N_fit .^ slope

p_scaling = plot(N_fit, R2_fit,
    label = "Fit: \$R^2_{ee} \\propto N^{$(round(slope,digits=3))},  \\nu=$(round(nu_fit,digits=3))\$",
    linestyle = :dash, color = :red;
    plot_kwargs()...,
    xscale = :log10, yscale = :log10,
    xlabel = L"N", ylabel = L"\langle R^2_{ee} \rangle",
    title  = "Scaling of end-to-end radius (T = 4.98)")

scatter!(p_scaling, N_vals_mc, R2_means,
    yerror = R2_errs,
    label  = "MC simulation",
    color  = :steelblue, markersize = 7)

savefig(p_scaling, joinpath(figdir, "scaling_R2ee.pdf"))


############################################

# Exercise 3
# idea: start at origin (0,0)
# then random aus ((+1,0), (-1,0), (0,+1), (0,-1))

function random_walk(N,n)

    arr = [[1,0], [-1,0], [0,1], [0,-1]]

    #config_list = []
    distance_list = []

    for _ in 1:n

        config = [0,0]

        for _ in 1:N 
            
            zufall = rand(arr)

            config = config .+ zufall

        end

        #push!(config_list,config)

        distance = config[1]^2 + config[2]^2
        push!( distance_list, distance)

    end 

    return sqrt(mean(distance_list))

end


N_vals = [200,400,600,800,1000]
n_val = 10_000

R_vals = [random_walk(N,n_val)[1] for N in N_vals]

p_random_walk = plot(N_vals, R_vals, title="Random walk running $n_val times", label = L"\bar{R}" , xaxis=:log, yaxis=:log, xlabel = L"N", ylabel = L"\bar{R}"; plot_kwargs()...);
plot!(p_random_walk,N_vals, N_vals .^(0.5), label = L"N^{0.5}", linewidth = 4,  linestyle=:dash );

savefig(p_random_walk, joinpath(figdir, "random_walk.pdf"))

############################################
#exercise 4

f(M,R,T,v) = 4*pi*(M/(2*pi*R*T))^(3/2) *v^2*exp(-M*v^2/(2*R*T))

# M = 60.08*10^(-3)kg/mol
# R = 8.314 H/(mol*K)
# T = 1450 K 

v = 0:0.1:2_000

plot_bolzmann = plot(v, f.(60.08*10^(-3), 8.314,1450,v), label = L"f(v)"; plot_kwargs()..., top_margin = 5 * Plots.mm)

savefig(plot_bolzmann, joinpath(figdir, "Bolzmann.pdf"))