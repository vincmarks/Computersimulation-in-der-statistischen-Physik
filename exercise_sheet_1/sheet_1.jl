# loading the packages
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


# Pkg.add("Distributions")
# Pkg.add("Plots")
# Pkg.add("LaTeXStrings")
# Pkg.add("Random")

using Random
using Distributions
using Plots
using LaTeXStrings

figdir = joinpath(@__DIR__, "Plots")
#######################
# Exercise 1
#######################
# (a)
x = rand(Uniform(0, 1), 100_000)

# (b)
function average(x)
    return sum(x) / length(x)
end

function variance(x)
    m = average(x)
    return sum((x .- m) .^ 2) / length(x) 
end

function error(x)
    return sqrt(variance(x) / length(x))
end

average(x), variance(x), error(x)

# (c)

fig_hist = histogram(x, bins=100, normalize=false, label="Data", title="Histogram of Uniform(0,1) Data", xlabel="Value", ylabel="Frequency")
savefig(fig_hist, joinpath(figdir, "histogram.pdf"))

# (d)
d_x = rand(Uniform(0, 1), 100_000)
d_y = rand(Uniform(0, 1), 100_000)
fig = scatter(d_x, d_y, label="Data", title="Scatter Plot of Uniform(0,1) Data", xlabel="X", ylabel="Y")
savefig(fig, joinpath(figdir, "scatter_plot.pdf"))

# (e)
function random_number_generator(seed, a, c, m, n)
    x = zeros(Int64, n)
    x[1] = seed
    for i in 2:n
        x[i] = (a * x[i-1] + c) % m #x_{n+1} = (a * x_n + c) mod m
    end
    return x/m # normalize to [0,1]
end

random_number_generator(1, 1664525, 1013904223, 2^32, 10)  # for example 10 random numbers

#######################
# Exercise 2
#######################
#(a)

function pi(n)
    total_points = 0

    for i in 1:n
        x = rand(Uniform(0, 1))
        y = rand(Uniform(0, 1))
        if x^2 + y^2 <= 1
            total_points += 1
        end
    end
    return (total_points / n) * 4
end


V = [10_000, 100_000, 1_000_000, 10_000_000]

for n in V
    println("Estimated Pi with $n points: ", pi(n))
end

function error_pi()
   
    num_trials = 100
    pi_true = π
    average_errors = Float64[]

    for n in V
        errors = Float64[]
        for trial in 1:num_trials
            pi_est = pi(n)
            error_val = abs(pi_est - pi_true)
            push!(errors, error_val)
        end
        avg_error = average(errors)
        push!(average_errors, avg_error)
        println("N=$n: Average Error = $avg_error")
    end

    # Plot N vs Average Error 
    fig_error = plot(V, average_errors, xscale=:log10, yscale=:log10, 
                label="Average Absolute Error", 
                 title="N vs. average absolute error",
                 xlabel="N",
                 ylabel="Average Absolute Error",
                marker=:circle, markersize=5, linewidth=2)

    return savefig(fig_error, joinpath(figdir, "error_plot.pdf"))

end

error_pi()


function visualization_pi_estimation(n)
    x_inside = Float64[]
    y_inside = Float64[]
    x_outside = Float64[]
    y_outside = Float64[]

    for i in 1:n
        x = rand()
        y = rand()
        if x^2 + y^2 <=1
            push!(x_inside,x)
            push!(y_inside,y)
        else
            push!(x_outside,x)
            push!(y_outside,y)
        end
    end

    scatter(x_inside, y_inside, label="Inside Circle", color=:blue, marker=:circle, markersize=3)
    scatter!(x_outside, y_outside, label="Outside Circle", color=:red, marker=:circle, markersize=3)
    x = range(0,1,100)
    y = sqrt.(1 .- x.^2)
    plot_pi_estimate =  scatter!(x, y, label=L"y = \sqrt{1 - x^2}", color=:black, linewidth=2)

    return savefig(plot_pi_estimate, joinpath(figdir, "pi_estimation_visualization.pdf"))
end

visualization_pi_estimation(100_000)