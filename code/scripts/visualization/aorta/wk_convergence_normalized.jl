# ==========================================================================================
# Plot normalized convergence of outlet pressures and outlet flows with varying particle spacings
# ==========================================================================================
# This script plots pressure and flow waveforms at each outlet for the last cycle,
# normalized by the finest resolution reference (dp = 0.00035).
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
reference_spacing = 0.00035
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
# ==== Plot normalized convergence
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

function get_last_cycle_interpolated(values, times_total, t_final, T, times_cycle)
    n_cycles = Int(floor(t_final / T))
    if n_cycles != 7
        @warn "Expected 7 cycles, got $n_cycles cycles. t_final: $t_final"
    end

    t_start = (n_cycles - 1) * T
    t_end = n_cycles * T

    cycle_mask = (times_total .>= t_start) .& (times_total .<= t_end)
    if sum(cycle_mask) == 0
        return nothing
    end

    t_cycle_data = times_total[cycle_mask] .- t_start
    cycle_data = values[cycle_mask]

    return [begin
                idx_interp = searchsortedfirst(t_cycle_data, t)
                if idx_interp > length(t_cycle_data)
                    cycle_data[end]
                elseif idx_interp == 1
                    cycle_data[1]
                else
                    t1, t2 = t_cycle_data[idx_interp - 1], t_cycle_data[idx_interp]
                    v1, v2 = cycle_data[idx_interp - 1], cycle_data[idx_interp]
                    v1 + (v2 - v1) * (t - t1) / (t2 - t1)
                end
            end
            for t in times_cycle]
end

function normalize_by_reference(values, ref_values; atol=1e-10)
    normalized = similar(values)
    for i in eachindex(values)
        denom = ref_values[i]
        normalized[i] = abs(denom) > atol ? values[i] / denom : NaN
    end
    return normalized
end

# function relative_error(values, ref_values; atol=1e-10)
#     eps = similar(values)
#     for i in eachindex(values)
#         denom = abs(ref_values[i])
#         eps[i] = denom > atol ? abs(values[i] - ref_values[i]) / denom : NaN
#     end
#     return eps
# end

# function trapz(t, y)
#     length(t) < 2 && return 0.0
#     integral = 0.0
#     for i in 1:(length(t) - 1)
#         dt = t[i + 1] - t[i]
#         integral += 0.5 * (y[i] + y[i + 1]) * dt
#     end
#     return integral
# end

# function cycle_average_error(t, eps)
#     mask = isfinite.(eps)
#     sum(mask) < 2 && return NaN
#     t_valid = t[mask]
#     eps_valid = eps[mask]
#     duration = t_valid[end] - t_valid[1]
#     duration <= 0 && return NaN
#     return trapz(t_valid, eps_valid) / duration
# end

# function normalized_l2_error(t, values, ref_values)
#     mask = isfinite.(values) .& isfinite.(ref_values)
#     sum(mask) < 2 && return NaN
#     t_valid = t[mask]
#     val_valid = values[mask]
#     ref_valid = ref_values[mask]

#     numerator = trapz(t_valid, (val_valid .- ref_valid) .^ 2)
#     denominator = trapz(t_valid, ref_valid .^ 2)

#     denominator <= 0 && return NaN
#     return sqrt(numerator / denominator)
# end

# function fit_loglog_slope(x, y)
#     mask = isfinite.(x) .& isfinite.(y) .& (x .> 0) .& (y .> 0)
#     sum(mask) < 2 && return NaN

#     lx = log10.(x[mask])
#     ly = log10.(y[mask])
#     lx_mean = mean(lx)
#     ly_mean = mean(ly)
#     denom = sum((lx .- lx_mean) .^ 2)
#     denom <= 0 && return NaN

#     return sum((lx .- lx_mean) .* (ly .- ly_mean)) / denom
# end

# Get outlet order from first param_sim - assume same for all
param_sim = results_all[1].param_sim
outlet_order = sort([key for (key, val) in param_sim.boundaries if key != "inflow"],
                    by=key -> param_sim.boundaries[key].id)

# Find reference dataset
reference_idx = findfirst(x -> isapprox(x.particle_spacing, reference_spacing; atol=1e-12),
                          results_all)
reference_idx === nothing &&
    error("Reference spacing $(reference_spacing) not found in particle_spacings.")
reference_data = results_all[reference_idx]

# Update theme with colormap-based palette for automatic color cycling
update_theme!(palette=(color=cgrad(:redblue, length(particle_spacings), rev=true,
                                   categorical=true),))

# ------------------------------------------------------------------------------------------
# ==== Normalized pressure plot
fig = Figure(size=(1000, 1200 / 5) .* 0.9)

axs = []
for (col, key) in enumerate(outlet_order)
    ax = Axis(fig[1, col],
              xlabel="",
              ylabel=col == 1 ? "p / p_ref [-]" : "",
              title=format_outlet_name(key),
              titlesize=12)
    push!(axs, ax)
end

for (col, key) in enumerate(outlet_order)
    outlet_id = param_sim.boundaries[key].id

    # Reference waveform for this outlet (finest resolution)
    p_ref = reference_data.results.p[outlet_id, :]
    p_ref_interp = get_last_cycle_interpolated(p_ref, reference_data.times_total,
                                               reference_data.t_final, T, times_cycle)

    if p_ref_interp === nothing
        @warn "No reference pressure data found for outlet $(key)."
        continue
    end

    for data in results_all
        p_outlet = data.results.p[outlet_id, :]
        p_interp = get_last_cycle_interpolated(p_outlet, data.times_total, data.t_final,
                                               T, times_cycle)
        p_interp === nothing && continue

        p_norm = normalize_by_reference(p_interp, p_ref_interp)

        lines!(axs[col], times_cycle, p_norm,
               label=format_dp_label(data.particle_spacing), linewidth=1.5)

        p_mean_ref = reference_data.results.p_mean[outlet_id, end]
        p_mean = data.results.p_mean[outlet_id, end]
        p_mean_norm = abs(p_mean_ref) > 1e-10 ? p_mean / p_mean_ref : NaN
        hlines!(axs[col], p_mean_norm, linestyle=:dash, linewidth=1.0)
    end

    xlims!(axs[col], low=0)
    ylims!(axs[col], low=0.7, high=1.1)
end

for i in 2:length(axs)
    hideydecorations!(axs[i], grid=false)
end

Legend(fig[1, length(outlet_order) + 1], axs[1], fontsize=10, halign=:left, valign=:center,
       tellheight=true)

Label(fig[2, :], "Time [s]", fontsize=14)
colgap!(fig.layout, 5)

resize_to_layout!(fig)

if save_fig
    dir = joinpath(fig_dir(), "aorta")
    mkpath(dir)
    save(joinpath(dir,
                  "wk_convergence_normalized_$(subject)_$(fsi ? "elastic" : "rigid").pdf"),
         fig)
else
    display(fig)
end
fig

# # ------------------------------------------------------------------------------------------
# # ==== Relative pressure error analysis (vs. reference solution)
# comparison_data = [d for d in results_all
#                    if !isapprox(d.particle_spacing, reference_spacing; atol=1e-12)]
# sort!(comparison_data, by=d -> d.particle_spacing)

# comparison_spacings = [d.particle_spacing for d in comparison_data]
# n_outlets = length(outlet_order)
# n_comp = length(comparison_data)

# eps_cycle_mean = fill(NaN, n_outlets, n_comp)
# eps_l2 = fill(NaN, n_outlets, n_comp)

# # 1) Relative error over time epsilon(t)
# fig_error_time = Figure(size=(1000, 1200 / 5) .* 0.9)
# axs_error_time = []

# for (col, key) in enumerate(outlet_order)
#     ax = Axis(fig_error_time[1, col],
#               xlabel="",
#               ylabel=col == 1 ? "ε(t) = |p - p_ref| / |p_ref| [-]" : "",
#               title=format_outlet_name(key),
#               titlesize=12)
#     push!(axs_error_time, ax)
# end

# for (col, key) in enumerate(outlet_order)
#     outlet_id = param_sim.boundaries[key].id

#     p_ref = reference_data.results.p[outlet_id, :]
#     p_ref_interp = get_last_cycle_interpolated(p_ref, reference_data.times_total,
#                                                reference_data.t_final, T, times_cycle)
#     p_ref_interp === nothing && continue

#     for (idx, data) in enumerate(comparison_data)
#         p_outlet = data.results.p[outlet_id, :]
#         p_interp = get_last_cycle_interpolated(p_outlet, data.times_total, data.t_final,
#                                                T, times_cycle)
#         p_interp === nothing && continue

#         eps_t = relative_error(p_interp, p_ref_interp)
#         lines!(axs_error_time[col], times_cycle, eps_t,
#                label=format_dp_label(data.particle_spacing), linewidth=1.5)

#         eps_cycle_mean[col, idx] = cycle_average_error(times_cycle, eps_t)
#         eps_l2[col, idx] = normalized_l2_error(times_cycle, p_interp, p_ref_interp)
#     end

#     xlims!(axs_error_time[col], low=0)
# end

# for i in 2:length(axs_error_time)
#     hideydecorations!(axs_error_time[i], grid=false)
# end

# Legend(fig_error_time[1, length(outlet_order) + 1], axs_error_time[1], fontsize=10,
#        halign=:left, valign=:center, tellheight=true)
# Label(fig_error_time[2, :], "Time [s]", fontsize=14)
# colgap!(fig_error_time.layout, 5)
# resize_to_layout!(fig_error_time)

# if save_fig
#     save(joinpath(fig_dir(), "aorta",
#                   "wk_convergence_error_time_$(subject)_$(fsi ? "elastic" : "rigid").pdf"),
#          fig_error_time)
# else
#     fig_error_time
# end

# # 2) Cycle-averaged relative error
# fig_error_cycle = Figure(size=(900, 500))
# ax_error_cycle = Axis(fig_error_cycle[1, 1],
#                       xlabel="Δx [m]",
#                       ylabel="ε_cycle [-]",
#                       title="Cycle-averaged relative pressure error",
#                       xscale=log10,
#                       yscale=log10)

# for (out_idx, key) in enumerate(outlet_order)
#     yvals = eps_cycle_mean[out_idx, :]
#     lines!(ax_error_cycle, comparison_spacings, yvals,
#            label=format_outlet_name(key), linewidth=1.4)
#     scatter!(ax_error_cycle, comparison_spacings, yvals, markersize=7)
# end

# axislegend(ax_error_cycle; position=:rb, labelsize=10)
# resize_to_layout!(fig_error_cycle)

# if save_fig
#     save(joinpath(fig_dir(), "aorta",
#                   "wk_convergence_error_cycle_$(subject)_$(fsi ? "elastic" : "rigid").pdf"),
#          fig_error_cycle)
# else
#     fig_error_cycle
# end

# # 3) Log-log convergence plot with normalized L2 error
# fig_l2 = Figure(size=(900, 500))
# ax_l2 = Axis(fig_l2[1, 1],
#              xlabel="Δx [m]",
#              ylabel="ε_L2 [-]",
#              title="Log-log convergence (normalized L2 pressure error)",
#              xscale=log10,
#              yscale=log10)

# for (out_idx, key) in enumerate(outlet_order)
#     yvals = eps_l2[out_idx, :]
#     lines!(ax_l2, comparison_spacings, yvals,
#            label=format_outlet_name(key), linewidth=1.2)
#     scatter!(ax_l2, comparison_spacings, yvals, markersize=7)
# end

# # Also show mean error across outlets as a thicker reference curve
# eps_l2_mean = [mean(filter(isfinite, eps_l2[:, i])) for i in 1:n_comp]
# lines!(ax_l2, comparison_spacings, eps_l2_mean,
#        color=:black, linewidth=2.5, label="Mean")
# scatter!(ax_l2, comparison_spacings, eps_l2_mean,
#          color=:black, markersize=9)

# slope = fit_loglog_slope(comparison_spacings, eps_l2_mean)
# if isfinite(slope)
#     text!(ax_l2, comparison_spacings[1], eps_l2_mean[1],
#           text="slope ≈ $(round(slope, digits=2))",
#           align=(:left, :bottom), fontsize=11, color=:black)
# end

# axislegend(ax_l2; position=:rb, labelsize=10)
# resize_to_layout!(fig_l2)

# if save_fig
#     save(joinpath(fig_dir(), "aorta",
#                   "wk_convergence_error_l2_loglog_$(subject)_$(fsi ? "elastic" : "rigid").pdf"),
#          fig_l2)
# else
#     fig_l2
# end
