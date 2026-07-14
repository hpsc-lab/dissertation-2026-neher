using SimulationSetup
using CSV, DataFrames
using CairoMakie

include("../theme.jl")

param1 = length(ARGS) >= 1 ? (ARGS[1]) : ("F09")

# Configuration
save_fig = true
subject = param1

# Hardcoded versions to compare
version1 = :normotensive
version2 = :hypertensive
version3 = :exercise

# ==========================================================================================
# ==== Load CSV data for both versions

# Version 1
csv_file1 = joinpath(data_dir(), "flow_rates", "volume_rate_$(version1)_$(subject).csv")
if !isfile(csv_file1)
    error("CSV file not found: $csv_file1\nPlease run volume_rate.jl for version $version1 first.")
end
data1 = CSV.read(csv_file1, DataFrame)

# Version 2
csv_file2 = joinpath(data_dir(), "flow_rates", "volume_rate_$(version2)_$(subject).csv")
if !isfile(csv_file2)
    error("CSV file not found: $csv_file2\nPlease run volume_rate.jl for version $version2 first.")
end
data2 = CSV.read(csv_file2, DataFrame)

# Version 3 (only for F09)
if subject == "F09"
    csv_file3 = joinpath(data_dir(), "flow_rates",
                         "volume_rate_$(version3)_$(subject).csv")
    if !isfile(csv_file3)
        error("CSV file not found: $csv_file3\nPlease run volume_rate.jl for version $version3 first.")
    end
    data3 = CSV.read(csv_file3, DataFrame)
end

# ==========================================================================================
# ==== Extract data and remove missing values

# Helper function to remove missing values
function remove_missing_values(times, values...)
    mask = .!ismissing.(times)
    times_clean = times[mask]
    values_clean = [v[mask] for v in values]
    return times_clean, values_clean...
end

# Version 1 data
times_rigid1, Q_analytic_rigid1, Q_total_out_rigid1,
volume_rigid_relative1 = remove_missing_values(data1.times_rigid,
                                               data1.Q_analytic_rigid,
                                               data1.Q_total_out_rigid,
                                               data1.volume_rigid_relative)
times_elastic1, Q_total_out_elastic1,
volume_elastic_relative1 = remove_missing_values(data1.times_elastic,
                                                 data1.Q_total_out_elastic,
                                                 data1.volume_elastic_relative)
Q_analytic_mean1 = data1.Q_analytic_mean[1]
Q_rigid_mean1 = data1.Q_rigid_mean[1]
Q_elastic_mean1 = data1.Q_elastic_mean[1]

# Version 2 data
times_rigid2, Q_analytic_rigid2, Q_total_out_rigid2,
volume_rigid_relative2 = remove_missing_values(data2.times_rigid,
                                               data2.Q_analytic_rigid,
                                               data2.Q_total_out_rigid,
                                               data2.volume_rigid_relative)
times_elastic2, Q_total_out_elastic2,
volume_elastic_relative2 = remove_missing_values(data2.times_elastic,
                                                 data2.Q_total_out_elastic,
                                                 data2.volume_elastic_relative)
Q_analytic_mean2 = data2.Q_analytic_mean[1]
Q_rigid_mean2 = data2.Q_rigid_mean[1]
Q_elastic_mean2 = data2.Q_elastic_mean[1]

# Version 3 data (only for F09)
if subject == "F09"
    times_rigid3, Q_analytic_rigid3, Q_total_out_rigid3,
    volume_rigid_relative3 = remove_missing_values(data3.times_rigid,
                                                   data3.Q_analytic_rigid,
                                                   data3.Q_total_out_rigid,
                                                   data3.volume_rigid_relative)
    times_elastic3, Q_total_out_elastic3,
    volume_elastic_relative3 = remove_missing_values(data3.times_elastic,
                                                     data3.Q_total_out_elastic,
                                                     data3.volume_elastic_relative)
    Q_analytic_mean3 = data3.Q_analytic_mean[1]
    Q_rigid_mean3 = data3.Q_rigid_mean[1]
    Q_elastic_mean3 = data3.Q_elastic_mean[1]
    t_final3 = max(maximum(times_rigid3), maximum(times_elastic3))
end

# Calculate time limits
t_final1 = max(maximum(times_rigid1), maximum(times_elastic1))
t_final2 = max(maximum(times_rigid2), maximum(times_elastic2))

# Plot colors
color_prescribed = Cycled(1)
color_rigid = Cycled(2)
color_elastic = Cycled(3)

# ==========================================================================================
# ==== Create comparison plot
set_theme!(my_thesis_theme)
fig = Figure(size=(1600, 800) .* 0.45)

# Left column: Version 1
ax1_1 = Axis(fig[1, 1], title="Normotensive")
ax1_2 = Axis(fig[2, 1])

# Version 1 - Plot 1: Volume flow rate
lines!(ax1_1, times_rigid1, Q_analytic_rigid1, label="q_in (prescribed)",
       color=color_prescribed, linewidth=2)
lines!(ax1_1, times_rigid1, Q_total_out_rigid1, label="∑q_out (rigid)",
       color=color_rigid, linewidth=2)
lines!(ax1_1, times_elastic1, Q_total_out_elastic1, label="∑q_out (elastic)",
       color=color_elastic, linewidth=2)

# Plot mean values
hlines!(ax1_1, Q_analytic_mean1, color=color_prescribed, linestyle=:dash, linewidth=1.5)
hlines!(ax1_1, Q_rigid_mean1, color=color_rigid, linestyle=:dash, linewidth=1.5)
hlines!(ax1_1, Q_elastic_mean1, color=color_elastic, linestyle=:dash, linewidth=1.5)

ax1_1.ylabel = "q [ml/s]"
ax1_1.yticks = LinearTicks(5)
ax1_1.xticks = 0:0.15:t_final1
hidexdecorations!(ax1_1, grid=false)
xlims!(ax1_1, low=0, high=t_final1)
ylims!(ax1_1, low=0, high=500)

# Version 1 - Plot 2: Volume change
lines!(ax1_2, times_rigid1, volume_rigid_relative1, label="rigid", color=color_rigid,
       linewidth=2)
lines!(ax1_2, times_elastic1, volume_elastic_relative1, label="elastic",
       color=color_elastic, linewidth=2)

ax1_2.xlabel = "Time [s]"
ax1_2.ylabel = "Volume change [%]"
ax1_2.yticks = LinearTicks(5)
ax1_2.xticks = 0:0.15:0.7
xlims!(ax1_2, low=0, high=t_final1)
ylims!(ax1_2, low=0, high=4.5)

# Middle column: Version 2
ax2_1 = Axis(fig[1, 2], title="Hypertensive")
ax2_2 = Axis(fig[2, 2])

# Version 2 - Plot 1: Volume flow rate
lines!(ax2_1, times_rigid2, Q_analytic_rigid2, label="Q_in (prescribed)",
       color=color_prescribed, linewidth=2)
lines!(ax2_1, times_rigid2, Q_total_out_rigid2, label="∑Q_out (rigid)",
       color=color_rigid, linewidth=2)
lines!(ax2_1, times_elastic2, Q_total_out_elastic2, label="∑Q_out (elastic)",
       color=color_elastic, linewidth=2)

# Plot mean values
hlines!(ax2_1, Q_analytic_mean2, color=color_prescribed, linestyle=:dash, linewidth=1.5)
hlines!(ax2_1, Q_rigid_mean2, color=color_rigid, linestyle=:dash, linewidth=1.5)
hlines!(ax2_1, Q_elastic_mean2, color=color_elastic, linestyle=:dash, linewidth=1.5)

ax2_1.ylabel = "q [ml/s]"
ax2_1.yticks = LinearTicks(5)
ax2_1.xticks = 0:0.15:t_final2
hidexdecorations!(ax2_1, grid=false)
hideydecorations!(ax2_1, grid=false)
xlims!(ax2_1, low=0, high=t_final2)
ylims!(ax2_1, low=0)
ylims!(ax2_1, low=0, high=500)

# Version 2 - Plot 2: Volume change
lines!(ax2_2, times_rigid2, volume_rigid_relative2, label="rigid", color=color_rigid,
       linewidth=2)
lines!(ax2_2, times_elastic2, volume_elastic_relative2, label="elastic",
       color=color_elastic, linewidth=2)

ax2_2.xlabel = "Time [s]"
ax2_2.ylabel = "Volume change [%]"
ax2_2.yticks = LinearTicks(5)
ax2_2.xticks = 0:0.15:t_final2
hideydecorations!(ax2_2, grid=false)
xlims!(ax2_2, low=0, high=t_final2)
ylims!(ax2_2, low=0, high=4.5)

# Right column: Version 3 (only for F09)
if subject == "F09"
    ax3_1 = Axis(fig[1, 3], title="Exercise")
    ax3_2 = Axis(fig[2, 3])

    # Version 3 - Plot 1: Volume flow rate
    lines!(ax3_1, times_rigid3, Q_analytic_rigid3, label="Q_in (prescribed)",
           color=color_prescribed, linewidth=2)
    lines!(ax3_1, times_rigid3, Q_total_out_rigid3, label="∑Q_out (rigid)",
           color=color_rigid, linewidth=2)
    lines!(ax3_1, times_elastic3, Q_total_out_elastic3, label="∑Q_out (elastic)",
           color=color_elastic, linewidth=2)

    # Plot mean values
    hlines!(ax3_1, Q_analytic_mean3, color=color_prescribed, linestyle=:dash, linewidth=1.5)
    hlines!(ax3_1, Q_rigid_mean3, color=color_rigid, linestyle=:dash, linewidth=1.5)
    hlines!(ax3_1, Q_elastic_mean3, color=color_elastic, linestyle=:dash, linewidth=1.5)

    ax3_1.ylabel = "q [ml/s]"
    ax3_1.yticks = LinearTicks(5)
    ax3_1.xticks = 0:0.15:t_final3
    hidexdecorations!(ax3_1, grid=false)
    # hideydecorations!(ax3_1, grid=false)
    xlims!(ax3_1, low=0, high=t_final3)
    ylims!(ax3_1, low=0, high=800)

    # Version 3 - Plot 2: Volume change
    lines!(ax3_2, times_rigid3, volume_rigid_relative3, label="rigid", color=color_rigid,
           linewidth=2)
    lines!(ax3_2, times_elastic3, volume_elastic_relative3, label="elastic",
           color=color_elastic, linewidth=2)

    ax3_2.xlabel = "Time [s]"
    ax3_2.ylabel = "Volume change [%]"
    ax3_2.yticks = LinearTicks(5)
    ax3_2.xticks = 0:0.15:t_final3
    # hideydecorations!(ax3_2, grid=false)
    xlims!(ax3_2, low=0, high=t_final3)
    ylims!(ax3_2, low=0, high=10)
end

# One legend per row, placed on the far right
legend_col = subject == "F09" ? 4 : 3
Legend(fig[1, legend_col], ax1_1)
Legend(fig[2, legend_col], ax1_2)

# Scale plot column widths with their time span so both columns use the same time-to-length mapping
colsize!(fig.layout, 1, Auto(t_final1))
colsize!(fig.layout, 2, Auto(t_final2))
if subject == "F09"
    colsize!(fig.layout, 3, Auto(t_final3))
    # colsize!(fig.layout, 4, Fixed(220))
end

colgap!(fig.layout, 10)
rowgap!(fig.layout, 12)

resize_to_layout!(fig)

# Optional: Save figure
if save_fig
    if subject == "F09"
        filename = "volume_rate_comparison_$(version1)_$(version2)_$(version3)_$(subject).pdf"
    else
        filename = "volume_rate_comparison_$(version1)_$(version2)_$(subject).pdf"
    end

    dir = joinpath(fig_dir(), "aorta")
    mkpath(dir)
    save(joinpath(dir, filename), fig)
end

fig
