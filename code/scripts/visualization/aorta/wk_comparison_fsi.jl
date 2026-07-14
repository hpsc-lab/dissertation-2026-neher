# ==========================================================================================
# Compare periodic state of outlet pressures across different versions
# ==========================================================================================
# This script compares the pressure waveforms at each outlet for different simulation versions:
# - v1.0.35: Normotensive
# - v1.0.36: Exercise
# - v1.0.37: Hypertensive
#
# One complete cardiac cycle is plotted for each version for comparison.
# ==========================================================================================

using SimulationSetup
using CSV, DataFrames
using CairoMakie

param1 = length(ARGS) >= 1 ? ARGS[1] : "F09"

include(pkgdir(SimulationSetup, "..", "scripts", "aorta", "velocity_functions.jl"))

set_config!(version=v"1.0.35")
initialize_code_version!()

save_fig = true
particle_spacing = 0.5e-3
subject = param1

# Define versions to compare
scenarios = [:normotensive, :hypertensive]
# version_labels = ["Normotensive", "Exercise", "Hypertensive"]
version_labels = ["Normotensive", "Hypertensive"]

# Define simulation parameters for each version
version_params = Dict(
    :normotensive => (T=0.75, p_syst=125.0, p_diast=75.0, stroke_volume_factor=1.0,
     v_peak_factor=1.0),
    :exercise => (T=0.4, p_syst=180.0, p_diast=85.0, stroke_volume_factor=1.2,
     v_peak_factor=2.2),
    :hypertensive => (T=0.55, p_syst=200.0, p_diast=120.0, stroke_volume_factor=1.05,
     v_peak_factor=1.3)
)

# ==========================================================================================
# ==== Load data for all versions
results_rigid = Dict()
results_fsi = Dict()

# Helper function to load data for a specific version
function load_version_data(scenario; fsi=false)

    # Build results directory path manually
    base_out = joinpath(out_dir(; result_variant=scenario))

    if fsi
        results_dir = joinpath(base_out, "aorta", "$subject", "elastic",
                               "dp_$(particle_spacing)_t_0.002", "full_cycle")
    else
        results_dir = joinpath(base_out, "aorta", "$subject", "rigid",
                               "dp_$(particle_spacing)", "full_cycle")
    end

    if !isdir(results_dir)
        @warn "Directory not found (skipping): $results_dir"
        return nothing
    end

    file = joinpath(results_dir, "resulting_pressures.csv")

    if !isfile(file)
        @warn "File not found: $file"
        return nothing
    end

    return CSV.read(file, DataFrame)
end

function extract_cycle_pressure(data)
    # Extract all pressure columns
    pressure_cols = [col for col in names(data) if startswith(col, "p_")]
    pressures = Dict()

    for col in pressure_cols
        pressures[col] = data[:, col]
    end

    return pressures
end

# Load data for each version
outlet_order = nothing
for (idx, scenario) in enumerate(scenarios)
    global outlet_order

    data_rigid = load_version_data(scenario; fsi=false)
    data_fsi = load_version_data(scenario; fsi=true)

    if data_rigid !== nothing
        results_rigid[scenario] = data_rigid
    end

    if data_fsi !== nothing
        results_fsi[scenario] = data_fsi
    end

    # Get outlet order from first scenario (all versions have same outlets)
    if outlet_order === nothing
        params = version_params[scenario]
        param_sim = SimulationParameters(subject; particle_spacing, T=params.T,
                                         q_prescribed=realistic_flow_ratios,
                                         stroke_volume_factor=params.stroke_volume_factor,
                                         v_peak_factor=params.v_peak_factor,
                                         p_syst=params.p_syst, p_diast=params.p_diast,
                                         L_eff=0.35)

        # Get outlet order from param_sim.boundaries, sorted by id
        # Filter out inflow (id = 0) and sort by id
        outlet_order = sort([key for (key, val) in param_sim.boundaries if key != "inflow"],
                            by=key -> param_sim.boundaries[key].id)
    end
end

# ==========================================================================================
# ==== Plot pressures for all versions
# ==========================================================================================
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

# Use outlet_order from boundaries
outlet_keys = outlet_order

# Convert mmHg to Pa (multiply by 133.322)
mmHg_to_Pa = 133.322

n_outlets = length(outlet_keys)

# Create figure with subplots (outlets in columns)
fig = Figure(size=(1000, 300) .* 0.9)

# Define colors for each version# Plot colors
colors = [Cycled(1), Cycled(3), Cycled(2)]

# Store line objects for legend
line_objects = []
line_labels = String[]

# Create grid of axes - one column per outlet
axs = []
for (col, key) in enumerate(outlet_keys)
    ax = Axis(fig[1, col],
              xlabel="",
              ylabel=col == 1 ? "Pressure [mmHg]" : "",
              title=format_outlet_name(key),
              titlesize=12,
              yticks=0:20:200)
    push!(axs, ax)
end

# Plot each outlet
for (col, key) in enumerate(outlet_keys)
    press_var = "p_outlet_$(key)_open_boundary_1"

    # Plot each version
    for (j, version) in enumerate(scenarios)
        if haskey(results_rigid, version)
            data = results_rigid[version]

            # Extract time and pressure
            time = data.time
            pressure_Pa = data[!, press_var]
            pressure_mmHg = pressure_Pa ./ mmHg_to_Pa

            # Plot
            l = lines!(axs[col], time, pressure_mmHg,
                       label="$(version_labels[j]) rigid",
                       color=colors[j],
                       linestyle=:dash)

            # Store line objects for legend (only once)
            if col == 1
                push!(line_objects, l)
                push!(line_labels, "$(version_labels[j]) rigid")
            end
        end

        if haskey(results_fsi, version)
            data = results_fsi[version]

            # Extract time and pressure
            time = data.time
            pressure_Pa = data[!, press_var]
            pressure_mmHg = pressure_Pa ./ mmHg_to_Pa

            # Plot
            l = lines!(axs[col], time, pressure_mmHg,
                       label="$(version_labels[j]) elastic",
                       color=colors[j],
                       linewidth=2)

            # Store line objects for legend (only once)
            if col == 1
                push!(line_objects, l)
                push!(line_labels, "$(version_labels[j]) elastic")
            end
        end
    end

    # Set y-axis limits
    ylims!(axs[col], low=20, high=160)
    xlims!(axs[col], low=0)
end

# Hide y-decorations for all but the first plot
for i in 2:length(axs)
    hideydecorations!(axs[i], grid=false)
end

# Add legend outside on the right
Legend(fig[1, length(outlet_keys) + 1], line_objects, line_labels,
       fontsize=10, halign=:left, valign=:center, tellheight=false)

# Add common x-label
Label(fig[2, :], "Time [s]", fontsize=14)
colgap!(fig.layout, 5)
resize_to_layout!(fig)

# Save figure
if save_fig
    dir = joinpath(fig_dir(), "aorta")
    mkpath(dir)
    save_path = joinpath(dir, "wk_comparison_pressures_fsi_$subject.pdf")
    save(save_path, fig)
    @info "Figure saved to: $save_path"
end

# display(fig)
