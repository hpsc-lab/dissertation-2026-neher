using SimulationSetup
using CSV, DataFrames
using CairoMakie, Statistics, LinearAlgebra

scenario = :normotensive

set_config!(version=v"1.0.35",scenario=scenario)
initialize_code_version!()

result_variant = scenario
save_fig = false

# Resolution
particle_spacing = 0.5e-3

# Patient identifier
subject = "F09"

include(pkgdir(SimulationSetup, "..", "scripts", "aorta", "velocity_functions.jl"))
include("../theme.jl")

if current_scenario() == :normotensive
    const T = 0.75 # period duration
elseif current_scenario() == :exercise # dynamic load
    const T = 0.4 # period duration
else # hypertensive peak
    const T = 0.55 # period duration
end
const omega = 2pi / T # Angular frequency

# ==========================================================================================
# ==== Load data
if current_scenario() == :normotensive
    p_syst = 125.0
    p_diast = 75.0
    stroke_volume_factor = 1.0
    v_peak_factor = 1.0
elseif current_scenario() == :exercise # dynamic load
    p_syst = 180.0
    p_diast = 85.0
    stroke_volume_factor = 1.2
    v_peak_factor = 2.2
else # hypertensive peak
    p_syst = 200.0
    p_diast = 120.0
    stroke_volume_factor = 1.05
    v_peak_factor = 1.3
end
param_sim = SimulationParameters(subject; particle_spacing, T,
                                 q_prescribed=realistic_flow_ratios,
                                 stroke_volume_factor=stroke_volume_factor,
                                 v_peak_factor=v_peak_factor,
                                 p_syst=p_syst, p_diast=p_diast, L_eff=0.35)

# Plot colors
color_prescribed = Cycled(1)
color_rigid = Cycled(2)
color_elastic = Cycled(3)

# ==========================================================================================
# ==== Load data from rigid and elastic simulations
results_dir_rigid = joinpath(out_dir(; result_variant), "aorta", "$subject", "rigid",
                             "dp_$(particle_spacing)", "full_cycle")
results_dir_elastic = joinpath(out_dir(; result_variant), "aorta", "$subject", "elastic",
                               "dp_$(particle_spacing)_t_0.002", "full_cycle")

file_rigid = joinpath(results_dir_rigid, "resulting_pressures.csv")
data_sim_rigid = CSV.read(file_rigid, DataFrame)

file_elastic = joinpath(results_dir_elastic, "resulting_pressures.csv")
data_sim_elastic = CSV.read(file_elastic, DataFrame)

times_rigid_ = data_sim_rigid[!, "time"]
times_elastic_ = data_sim_elastic[!, "time"]
t_final_rigid = last(times_rigid_)
t_final_elastic = last(times_elastic_)
times_rigid = times_rigid_[times_rigid_ .<= t_final_rigid]
times_elastic = times_elastic_[times_elastic_ .<= t_final_elastic]
t_final = max(t_final_rigid, t_final_elastic)

# ==========================================================================================
# ==== Preprocess simulation results
flow_rate_correction_factor = current_version() <= v"1.0.21" ? 1.3 : 1.0
Q_analytic_rigid = param_sim.subject_parameters.v_peak * param_sim.subject_parameters.A_in *
                   flow_rate_correction_factor * m3_to_ml() *
                   velocity_inlet_fourier.(times_rigid)
Q_analytic_elastic = param_sim.subject_parameters.v_peak *
                     param_sim.subject_parameters.A_in * flow_rate_correction_factor *
                     m3_to_ml() * velocity_inlet_fourier.(times_elastic)

# Process rigid results - sum of outlets
Q_outlets_rigid = zeros(param_sim.n_outlets, length(times_rigid))
for key in keys(param_sim.boundaries)
    key == "inflow" && continue
    Q = data_sim_rigid[!, "Q_outlet_$(key)_open_boundary_1"][1:length(times_rigid)]
    Q_outlets_rigid[param_sim.boundaries[key].id, :] .= Q
end
Q_total_out_rigid = vec(sum(Q_outlets_rigid, dims=1)) .* m3_to_ml()

# Process elastic results - sum of outlets
Q_outlets_elastic = zeros(param_sim.n_outlets, length(times_elastic))
for key in keys(param_sim.boundaries)
    key == "inflow" && continue
    Q = data_sim_elastic[!, "Q_outlet_$(key)_open_boundary_1"][1:length(times_elastic)]
    Q_outlets_elastic[param_sim.boundaries[key].id, :] .= Q
end
Q_total_out_elastic = vec(sum(Q_outlets_elastic, dims=1)) .* m3_to_ml()

# Calculate mean values
Q_analytic_mean = mean(Q_analytic_rigid)
Q_rigid_mean = mean(Q_total_out_rigid)
Q_elastic_mean = mean(Q_total_out_elastic)

# Extract total volume data and normalize to initial volume
volume_rigid = data_sim_rigid.total_volume_fluid_1[1:length(times_rigid)]
volume_elastic = data_sim_elastic.total_volume_fluid_1[1:length(times_elastic)]

# Calculate relative volume change in percent
volume_rigid_initial = minimum(volume_rigid)
volume_elastic_initial = minimum(volume_elastic)
volume_rigid_relative = (volume_rigid .- volume_rigid_initial) ./ volume_rigid_initial .*
                        100.0
volume_elastic_relative = (volume_elastic .- volume_elastic_initial) ./
                          volume_elastic_initial .* 100.0

# ==========================================================================================
# ==== Plot
set_theme!(my_thesis_theme)
fig = Figure(size=(1400, 800) .* 0.5)
g = fig[1, 1] = GridLayout()
ax1 = Axis(g[1, 1])
ax2 = Axis(g[2, 1])

# Plot 1: prescribed volume flow and sum of outlet flows
l1 = lines!(ax1, times_rigid, Q_analytic_rigid, label="q_in (prescribed)",
            color=color_prescribed, linewidth=2)
l2 = lines!(ax1, times_rigid, Q_total_out_rigid, label="∑q_out (rigid)",
            color=color_rigid, linewidth=2)
l3 = lines!(ax1, times_elastic, Q_total_out_elastic, label="∑q_out (elastic)",
            color=color_elastic, linewidth=2)

# Plot mean values as horizontal dashed lines
hlines!(ax1, Q_analytic_mean, color=color_prescribed, linestyle=:dash, linewidth=1.5)
hlines!(ax1, Q_rigid_mean, color=color_rigid, linestyle=:dash, linewidth=1.5)
hlines!(ax1, Q_elastic_mean, color=color_elastic, linestyle=:dash, linewidth=1.5)

Legend(g[1, 2], ax1)
ax1.ylabel = "Volume flow rate [ml/s]"
ax1.yticks = LinearTicks(5)
ax1.xticks = 0:0.15:t_final
hidexdecorations!(ax1, grid=false)

xlims!(ax1, low=0, high=t_final)
ylims!(ax1, low=0)

# Plot 2: relative volume change
l4 = lines!(ax2, times_rigid, volume_rigid_relative, label="rigid", color=color_rigid,
            linewidth=2)
l5 = lines!(ax2, times_elastic, volume_elastic_relative, label="elastic",
            color=color_elastic, linewidth=2)

Legend(g[2, 2], ax2)
ax2.xlabel = "Time [s]"
ax2.ylabel = "Volume change [%]"
ax2.yticks = LinearTicks(5)
ax2.xticks = 0:0.15:t_final

xlims!(ax2, low=0, high=t_final)
# ylims!(ax2, low=0)

rowgap!(g, 10)
resize_to_layout!(fig)

# Optional: Save figure
if save_fig
    save(joinpath(fig_dir(), "aorta", "volume_rate_$(current_scenario())_$(subject).pdf"),
         fig)
end

# ==========================================================================================
# ==== Save plot data to CSV
csv_filename = "volume_rate_$(current_scenario())_$(subject).csv"
csv_path = joinpath(data_dir(), "flow_rates")
mkpath(csv_path)

# Create DataFrame with all plot data
# Since times_rigid and times_elastic can have different lengths, we need to handle this
max_length = max(length(times_rigid), length(times_elastic))

# Helper function to pad arrays to same length with missing values
function pad_to_length(arr, target_length)
    if length(arr) == target_length
        return arr
    else
        return vcat(arr, fill(missing, target_length - length(arr)))
    end
end

# Create DataFrame with padded arrays
plot_data = DataFrame(times_rigid=pad_to_length(times_rigid, max_length),
                      Q_analytic_rigid=pad_to_length(Q_analytic_rigid, max_length),
                      Q_total_out_rigid=pad_to_length(Q_total_out_rigid, max_length),
                      volume_rigid_relative=pad_to_length(volume_rigid_relative,
                                                          max_length),
                      times_elastic=pad_to_length(times_elastic, max_length),
                      Q_total_out_elastic=pad_to_length(Q_total_out_elastic, max_length),
                      volume_elastic_relative=pad_to_length(volume_elastic_relative,
                                                            max_length),
                      Q_analytic_mean=fill(Q_analytic_mean, max_length),
                      Q_rigid_mean=fill(Q_rigid_mean, max_length),
                      Q_elastic_mean=fill(Q_elastic_mean, max_length))

CSV.write(joinpath(csv_path, csv_filename), plot_data)
println("Plot data saved to: $csv_path")

fig
