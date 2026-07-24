# ==========================================================================================
# Plot convergence of outlet pressures and outlet flows with varying particle spacings
# ==========================================================================================
# This script plots the pressure and flow waveforms at each outlet for the last cycle
# over different particle spacings to assess grid convergence.
#
# Data structure:
# - results.p: Matrix where rows are outlets and columns are time points
# - results.Q: Matrix where rows are outlets and columns are time points
# - times_total: Time vector covering all cycles
# - times_cycle: Time vector for one cycle (0 to T)
# ==========================================================================================

using SimulationSetup
using CSV, DataFrames
using Statistics
using CairoMakie

include(pkgdir(SimulationSetup, "..", "scripts", "aorta", "velocity_functions.jl"))

set_config!(version=v"1.0.35", scenario=:normotensive)
initialize_code_version!()

result_variant = current_scenario()

fsi = false
save_fig = true

particle_spacings = [0.001, 0.00075, 0.0005, 0.00035]
subject = "F09"

# ==========================================================================================
# ==== Load data
if current_scenario() == :normotensive
    const T = 0.75 # period duration
    p_syst = 125.0
    p_diast = 75.0
    stroke_volume_factor = 1.0
    v_peak_factor = 1.0
elseif current_scenario() == :exercise
    const T = 0.4 # period duration
    p_syst = 180.0
    p_diast = 85.0
    stroke_volume_factor = 1.2
    v_peak_factor = 2.2
else # hypertensive peak
    const T = 0.55 # period duration
    p_syst = 200.0
    p_diast = 120.0
    stroke_volume_factor = 1.05
    v_peak_factor = 1.3
end
const omega = 2pi / T # Angular frequency

# Load data for each particle spacing
results_all = []
times_cycle = collect(0.0:0.01:T)

for particle_spacing in particle_spacings
    param_sim = SimulationParameters(subject; particle_spacing, T,
                                     q_prescribed=realistic_flow_ratios,
                                     stroke_volume_factor=stroke_volume_factor,
                                     v_peak_factor=v_peak_factor,
                                     p_syst=p_syst, p_diast=p_diast, L_eff=0.35)

    if fsi
        results_dir = joinpath(out_dir(; result_variant), "aorta", "$subject", "elastic",
                               "dp_$(particle_spacing)_t_0.002")
    else
        results_dir = joinpath(out_dir(; result_variant), "aorta", "$subject", "rigid",
                               "dp_$(particle_spacing)")
    end
    file = joinpath(results_dir, "resulting_pressures.csv")
    data_sim = CSV.read(file, DataFrame)

    t_final = last(data_sim[!, "time"])
    times_total = data_sim[!, "time"][data_sim[!, "time"] .<= t_final]

    # Prepare simulation results for plotting
    results = prepare_plot_data(data_sim, param_sim, times_total; T=T, t_final=t_final)

    # Store results with metadata
    push!(results_all,
          (results=results, param_sim=param_sim, times_total=times_total,
           t_final=t_final, particle_spacing=particle_spacing))
end

function Q_prescr(t)
    flow_rate_correction_factor = 1.0
    param_sim = results_all[1].param_sim
    return param_sim.subject_parameters.v_peak * flow_rate_correction_factor *
           param_sim.subject_parameters.A_in * m3_to_ml() * velocity_inlet_fourier(t)
end

# ==========================================================================================
# ==== Plot convergence
include("../theme.jl")
set_theme!(my_thesis_theme)

# Function to format outlet names for display
function format_outlet_name(name)
    replacements = Dict(
        "left_subclavian" => "LSA",
        "right_common" => "RCCA",
        "left_common" => "LCCA",
        "right_subclavian" => "RSA",
        "thoracic" => "TA",
        "brachiocephalic" => "BCT"
    )
    return get(replacements, name, titlecase(replace(name, "_" => " ")))
end

# Function to format particle spacing for legend
function format_dp_label(dp)
    return "Δx = $((dp * 1e3))"
end

# Get outlet order from first param_sim - assume same for all
param_sim = results_all[1].param_sim
outlet_order = sort([key for (key, val) in param_sim.boundaries if key != "inflow"],
                    by=key -> param_sim.boundaries[key].id)

# Update theme with colormap-based palette for automatic color cycling
# update_theme!(palette=(color=cgrad(:darktest, length(particle_spacings),# rev=false,
#                                    categorical=true),))
update_theme!(palette=(color=cgrad(:redblue, length(particle_spacings), rev=true,
                                   categorical=true),))

# Create figure with subplots for each outlet
fig = Figure(size=(1000, 1200 / 5) .* 0.9)

# Create grid of axes - one column per outlet
axs = []
for (col, key) in enumerate(outlet_order)
    ax = Axis(fig[1, col],
              xlabel="",
              ylabel=col == 1 ? "Pressure [mmHg]" : "",
              title=format_outlet_name(key),
              titlesize=12)
    push!(axs, ax)
end

# Plot each outlet
for (col, key) in enumerate(outlet_order)
    outlet_id = param_sim.boundaries[key].id

    # Plot each particle spacing
    for (idx, data) in enumerate(results_all)
        results = data.results
        times_total = data.times_total
        t_final = data.t_final
        particle_spacing = data.particle_spacing

        # Calculate number of cycles
        n_cycles = Int(floor(t_final / T))
        if n_cycles != 7
            @warn "Expected 7 cycles, got $n_cycles cycles for particle spacing $particle_spacing. t_final: $t_final"
        end

        # Extract last cycle
        t_start = (n_cycles - 1) * T
        t_end = n_cycles * T

        # Get pressure data for this outlet
        p_outlet = results.p[outlet_id, :]

        # Find indices in times_total that correspond to the last cycle
        cycle_mask = (times_total .>= t_start) .& (times_total .<= t_end)

        if sum(cycle_mask) > 0
            t_cycle_data = times_total[cycle_mask] .- t_start
            p_cycle_data = p_outlet[cycle_mask]

            # Interpolate to times_cycle
            p_interpolated = [begin
                                  idx_interp = searchsortedfirst(t_cycle_data, t)
                                  if idx_interp > length(t_cycle_data)
                                      p_cycle_data[end]
                                  elseif idx_interp == 1
                                      p_cycle_data[1]
                                  else
                                      # Linear interpolation
                                      t1,
                                      t2 = t_cycle_data[idx_interp - 1],
                                           t_cycle_data[idx_interp]
                                      p1,
                                      p2 = p_cycle_data[idx_interp - 1],
                                           p_cycle_data[idx_interp]
                                      p1 + (p2 - p1) * (t - t1) / (t2 - t1)
                                  end
                              end
                              for t in times_cycle]

            # Plot with automatic color from palette
            lines!(axs[col], times_cycle, p_interpolated,
                   label=format_dp_label(particle_spacing), linewidth=1.5)
        end

        # Plot mean pressure as dashed horizontal line (from finest resolution)
        p_mean_outlet = data.results.p_mean[outlet_id, end]
        hlines!(axs[col], p_mean_outlet, linestyle=:dash, linewidth=1.0)
    end

    # # Plot prescribed mean pressure as dashed horizontal line
    # p_mean = (p_syst + 2 * p_diast) / 3
    # hlines!(axs[col], p_mean, color=:black, linestyle=:dash, linewidth=1.5)

    # Set x and y axis limits and ticks
    xlims!(axs[col], low=0)
    ylims!(axs[col], low=30, high=110)
    # ylims!(axs[col], low=0)
end

# Hide y-decorations for all but the first plot
for i in 2:length(axs)
    hideydecorations!(axs[i], grid=false)
end

# Add legend outside on the right
Legend(fig[1, length(outlet_order) + 1], axs[1], fontsize=10, halign=:left, valign=:center,
       tellheight=true)

# Add common x-label
Label(fig[2, :], "Time [s]", fontsize=14)
colgap!(fig.layout, 5)

resize_to_layout!(fig)

# Optional: Save figure
if save_fig
    dir = joinpath(fig_dir(), "aorta")
    mkpath(dir)
    save(joinpath(dir, "wk_convergence_$(subject)_$(fsi ? "elastic" : "rigid").pdf"), fig)
else
    fig
end

# ==========================================================================================
# ==== Plot flow convergence
fig_flow = Figure(size=(1000, 1200 / param_sim.n_outlets) .* 0.9)

# Create grid of axes - one column per outlet
axs_flow = []
for (col, key) in enumerate(outlet_order)
    ax = Axis(fig_flow[1, col],
              xlabel="",
              ylabel=col == 1 ? "Flow [ml/s]" : "",
              title=format_outlet_name(key),
              titlesize=12)
    push!(axs_flow, ax)
end

# Plot each outlet
for (col, key) in enumerate(outlet_order)
    outlet_id = param_sim.boundaries[key].id

    # Plot each particle spacing
    for data in results_all
        results = data.results
        times_total = data.times_total
        t_final = data.t_final
        particle_spacing = data.particle_spacing

        # Calculate number of cycles
        n_cycles = Int(floor(t_final / T))
        if n_cycles != 7
            @warn "Expected 7 cycles, got $n_cycles cycles for particle spacing $particle_spacing. t_final: $t_final"
        end

        # Extract last cycle
        t_start = (n_cycles - 1) * T
        t_end = n_cycles * T

        # Get flow data for this outlet
        Q_outlet = results.Q[outlet_id, :]

        # Find indices in times_total that correspond to the last cycle
        cycle_mask = (times_total .>= t_start) .& (times_total .<= t_end)

        if sum(cycle_mask) > 0
            t_cycle_data = times_total[cycle_mask] .- t_start
            Q_cycle_data = Q_outlet[cycle_mask]

            # Interpolate to times_cycle
            Q_interpolated = [begin
                                  idx_interp = searchsortedfirst(t_cycle_data, t)
                                  if idx_interp > length(t_cycle_data)
                                      Q_cycle_data[end]
                                  elseif idx_interp == 1
                                      Q_cycle_data[1]
                                  else
                                      # Linear interpolation
                                      t1,
                                      t2 = t_cycle_data[idx_interp - 1],
                                           t_cycle_data[idx_interp]
                                      Q1,
                                      Q2 = Q_cycle_data[idx_interp - 1],
                                           Q_cycle_data[idx_interp]
                                      Q1 + (Q2 - Q1) * (t - t1) / (t2 - t1)
                                  end
                              end
                              for t in times_cycle]

            # Plot with automatic color from palette
            lines!(axs_flow[col], times_cycle, Q_interpolated,
                   label=format_dp_label(particle_spacing), linewidth=1.5)
        end

        # Plot mean flow as dashed horizontal line
        Q_mean_outlet = data.results.Q_mean[outlet_id, end]
        hlines!(axs_flow[col], Q_mean_outlet, linestyle=:dash, linewidth=1.0)
    end

    # Set x and y axis limits
    xlims!(axs_flow[col], low=0)
end

# Hide y-decorations for all but the first and second plot
# for i in 3:length(axs_flow)
#     hideydecorations!(axs_flow[i], grid=false)
# end

# Add legend outside on the right
Legend(fig_flow[1, length(outlet_order) + 1], axs_flow[1], fontsize=10, halign=:left,
       valign=:center, tellheight=true)

# Add common x-label
Label(fig_flow[2, :], "Time [s]", fontsize=14)
colgap!(fig_flow.layout, 5)

resize_to_layout!(fig_flow)

# Optional: Save figure
if save_fig
    dir = joinpath(fig_dir(), "aorta")
    mkpath(dir)
    save(joinpath(dir, "wk_convergence_flow_$(subject)_$(fsi ? "elastic" : "rigid").pdf"),
         fig_flow)
else
    fig_flow
end
