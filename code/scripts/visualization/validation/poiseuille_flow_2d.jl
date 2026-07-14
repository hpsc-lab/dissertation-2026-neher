using SimulationSetup
using CairoMakie
using CSV, DataFrames

include(pkgdir(SimulationSetup, "..", "scripts", "plots", "theme.jl"))
set_theme!(my_thesis_theme)

save_fig = true

wall_distance = 0.001 # distance between top and bottom wall
flow_length = 0.004   # distance between inflow and outflow

fluid_density = 1000.0
reynolds_number = 50
pressure_drop = 0.1
dynamic_viscosity = sqrt(fluid_density * wall_distance^3 * pressure_drop /
                         (8 * flow_length * reynolds_number))

v_max = wall_distance^2 * pressure_drop / (8 * dynamic_viscosity * flow_length)

# Analytical velocity evolution given in eq. 16 (Zhang et al., 2025)
function poiseuille_velocity(y, t)

    # Base profile (stationary part)
    base_profile = (pressure_drop / (2 * dynamic_viscosity * flow_length)) * y *
                   (y - wall_distance)

    # Transient terms (Fourier series)
    transient_sum = 0.0

    for n in 0:10  # Limit to 10 terms for convergence
        coefficient = (4 * pressure_drop * wall_distance^2) /
                      (dynamic_viscosity * flow_length * pi^3 * (2 * n + 1)^3)

        sine_term = sin(pi * y * (2 * n + 1) / wall_distance)

        exp_term = exp(-((2 * n + 1)^2 * pi^2 * dynamic_viscosity * t) /
                       (fluid_density * wall_distance^2))

        transient_sum += coefficient * sine_term * exp_term
    end

    # Total velocity
    v_x = base_profile + transient_sum

    return v_x
end

# Load results
output_directory = joinpath(out_dir(), "validation", "open_boundaries",
                            "poiseuille_flow_2d")
data = CSV.read(joinpath(output_directory, "result_vx.csv"), DataFrame)

times = data[!, "time"]
times_ref = [0.1, 0.3, 0.6, 0.9, 2.0]
positions = range(0, wall_distance, length=100)
data_range = 2:98
data_indices = findall(t -> t in times_ref, times)
v_x_vector = [eval(Meta.parse(str)) for str in data[!, "v_x_fluid_1"]][data_indices]

# Calculate RMSEP error (eq. 17, Zhang et al., 2025)
rmsep_run = Float64[]
for (i, t) in enumerate(times_ref)
    N = length(data_range)
    res = sum(data_range, init=0) do j
        v_x = v_x_vector[i][j]

        v_analytical = -poiseuille_velocity(positions[j], t)

        # Avoid dividing by zero
        v_analytical < sqrt(eps()) && return 0.0

        rel_err = (v_analytical - v_x) / v_analytical

        return rel_err^2 / N
    end

    push!(rmsep_run, sqrt(res) * 100)
end

# RMSEP error (%) received by Zhang et al. (2025)
rmsep_reference = [1.81, 0.95, 0.67, 0.86, 1.22]

# recieved from my validation
rmsep_run = [0.585675, 0.481812, 0.785568, 1.12716, 1.62536]

# First plot: RMSEP error comparison
fig1 = Figure(size=(800, 600) .* 0.9)
ax1 = Axis(fig1[1, 1],
           xlabel="t", ylabel="RMSEP error (%)",
           limits=(0, 2.05, 0, 4))
scatter!(ax1, times_ref, rmsep_run, markersize=20, label="TrixiP")
scatter!(ax1, times_ref, rmsep_reference, marker=:x, markersize=20,
         label="Zhang et al. (2025)")
axislegend(ax1)

plot_range = range(0, wall_distance, length=50)
v_x_plot = view(stack(v_x_vector), 1:2:100, :)
label_ = "numerical (" .* ["0.1" "0.3" "0.6" "0.9" "2.0"] .* " s)"

# Second plot: Velocity profiles
fig2 = Figure(size=(850, 500) .* 0.8)
ax2 = Axis(fig2[1, 1],
           xlabel="y position (m)", ylabel="x velocity (m/s)",
           limits=(0, wall_distance, -0.002, 0.014))

# Plot simulation results
for i in 1:length(times_ref)
    scatter!(ax2, plot_range, v_x_plot[:, i], markersize=10, label=label_[i], marker=:x)
end

# Plot analytical solutions
for t in times_ref
    label__ = t == 2.0 ? "analytical" : nothing
    lines!(ax2, plot_range, (y) -> -poiseuille_velocity(y, t),
           linewidth=2, linestyle=:dash, color=:black, label=label__)
end

Legend(fig2[1, 2], ax2; tellwidth=true)

if save_fig
    dir = joinpath(fig_dir(), "validation")
    mkpath(dir)
    save(joinpath(dir, "poiseuille_flow.pdf"), fig2)
else
    fig2
end
