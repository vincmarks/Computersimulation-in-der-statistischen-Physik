 # wsl, vmd
 
# loading the packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


Pkg.add("Distributions")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
Pkg.add("Random")

using Random
using Distributions
using Plots
using LaTeXStrings


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

