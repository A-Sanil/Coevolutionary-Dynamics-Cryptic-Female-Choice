# Julia script to create visualization images for the three offspring functions
# Shows how each function should look: parabola, linear (monotonic), and asymptotic

using Plots

# Define the functions matching the Julia code
function offspring_parabola(score::Float64)
    peak = 4.0
    opt = 0.65
    width = 0.70
    scale = 1 - ((score - opt) / width)^2
    base = peak * max(0.0, scale)
    return max(1.0, base)
end

function offspring_asymptote(score::Float64)
    max_children = 4.0
    k = 4.2
    return max_children * (1 - exp(-k * score))
end

function offspring_linear(score::Float64)
    min_offspring = 1.0
    max_offspring = 4.0
    return min_offspring + (max_offspring - min_offspring) * score
end

# Create data for plotting
scores = 0.0:0.001:1.0
parabola_vals = [offspring_parabola(s) for s in scores]
asymptote_vals = [offspring_asymptote(s) for s in scores]
linear_vals = [offspring_linear(s) for s in scores]

# Create combined plot
p_combined = plot(scores, parabola_vals, 
    label="Parabola", 
    linewidth=2, 
    color=:orange,
    xlabel="Compatibility Score",
    ylabel="Expected Offspring Count",
    title="Offspring Functions: Expected Shapes",
    ylims=(0, 5),
    legend=:topright,
    dpi=300,
    size=(1200, 800))

plot!(p_combined, scores, asymptote_vals, 
    label="Asymptotic", 
    linewidth=2, 
    linestyle=:dash,
    color=:green)

plot!(p_combined, scores, linear_vals, 
    label="Linear (Monotonic)", 
    linewidth=2, 
    linestyle=:dot,
    color=:blue)

savefig(p_combined, "function_shapes_visualization.png")
println("Saved: function_shapes_visualization.png")

# Create individual plots
# Parabola
p1 = plot(scores, parabola_vals, 
    linewidth=2, 
    color=:orange,
    xlabel="Compatibility Score",
    ylabel="Expected Offspring Count",
    title="Parabola Function",
    ylims=(0, 5),
    dpi=300,
    size=(800, 600))
savefig(p1, "function_parabola.png")
println("Saved: function_parabola.png")

# Asymptotic
p2 = plot(scores, asymptote_vals, 
    linewidth=2, 
    linestyle=:dash,
    color=:green,
    xlabel="Compatibility Score",
    ylabel="Expected Offspring Count",
    title="Asymptotic Function",
    ylims=(0, 5),
    dpi=300,
    size=(800, 600))
savefig(p2, "function_asymptotic.png")
println("Saved: function_asymptotic.png")

# Linear (Monotonic)
p3 = plot(scores, linear_vals, 
    linewidth=2, 
    linestyle=:dot,
    color=:blue,
    xlabel="Compatibility Score",
    ylabel="Expected Offspring Count",
    title="Linear (Monotonically Increasing) Function",
    ylims=(0, 5),
    dpi=300,
    size=(800, 600))
savefig(p3, "function_linear.png")
println("Saved: function_linear.png")

println("All function visualization images created successfully!")

