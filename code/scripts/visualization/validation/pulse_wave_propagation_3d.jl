using SimulationSetup
using CairoMakie
using CSV, DataFrames

include(pkgdir(SimulationSetup, "..", "scripts", "plots", "theme.jl"))
set_theme!(my_thesis_theme)

save_fig = true
sound_speed_factor = 50
# Load results
output_directory = joinpath(out_dir(), "validation", "fsi",
                            "pulse_wave_propagation_3d_c_$sound_speed_factor")

# ======================================================================================
# ==== Read results
data = CSV.read(joinpath(output_directory, "values.csv"), DataFrame)

times = data[!, "time"]
times_ref = [1.6, 3.2, 4.8, 6.4, 8.0] .* 1e-3
data_indices = findall(t -> t in times_ref, round.(times, digits=4))
radial_displacements = [eval(Meta.parse(str))
                        for str in data[!, "radial_displacement_structure_1"]][data_indices]
pressures_along_axis = [eval(Meta.parse(str))
                        for str in data[!, "pressure_along_axis_fluid_1"]][data_indices]

# ======================================================================================
# ==== Read reference results

# corresponds to 1.6ms, 3.2ms, 4.8ms, 6.4ms, 8ms
curves = ["grey", "red", "blue", "green", "purple"]

radial_displacements_ref = Vector{DataFrame}(undef, length(curves))
pressures_along_axis_ref = Vector{DataFrame}(undef, length(curves))

for (i, curve) in enumerate(curves)
    radial_displacements_ref[i] = CSV.read(joinpath(data_dir(), "reference_data", "fsi",
                                                    "displacement",
                                                    "curve_$curve.csv"), DataFrame)
    pressures_along_axis_ref[i] = CSV.read(joinpath(data_dir(), "reference_data", "fsi",
                                                    "pressure",
                                                    "curve_$curve.csv"), DataFrame)
end

# ======================================================================================
# ==== Plot

plot_range = range(0, 0.1, length=length(first(radial_displacements)))
displacements = view(stack(radial_displacements), :, :)
pressures = view(stack(pressures_along_axis), :, :)
label_ = ["1.6", "3.2", "4.8", "6.4", "8"] .* " ms"

# Create figure with two subplots
fig = Figure(size=(900, 350))

# First subplot - Radial displacement
ax1 = Axis(fig[1, 1],
           xlabel="Distance (m)",
           ylabel="Radial displacement (m)",
           limits=(0, 0.1, -1e-4, 5e-4))

for i in 1:length(times_ref)
    lines!(ax1, plot_range, displacements[:, i],
           label=label_[i],
           linewidth=3)

    # reference scatter
    scatter!(ax1,
             radial_displacements_ref[i][!, 1] ./ 100,
             radial_displacements_ref[i][!, 2] ./ 100,
             markersize=6, marker=:x,
             label="Ref " * label_[i])
end

# Second subplot - Pressure
ax2 = Axis(fig[1, 2],
           xlabel="Distance (m)",
           ylabel="Pressure along the centerline",
           limits=(0, 0.1, -1500, 7500))

for i in 1:length(times_ref)
    lines!(ax2, plot_range, pressures[:, i],
           label=label_[i],
           linewidth=3)

    # reference scatter
    scatter!(ax2,
             pressures_along_axis_ref[i][!, 1] ./ 100,
             pressures_along_axis_ref[i][!, 2] .* 1000,
             markersize=6, marker=:x,
             label="Ref " * label_[i])
end

# Add legend to the right
Legend(fig[1, 3], ax2, framevisible=true)

if save_fig
    dir = joinpath(fig_dir(), "validation")
    mkpath(dir)
    save(joinpath(dir, "flexible_pipe_validation_$sound_speed_factor.pdf"), fig)
else
    fig
end
