# Plot validation of Windkessel pressure models (Fig. 2.6)
# Includes validation/pressure_model/windkessel_model.jl to generate sol_2_element and sol_3_element
# Compares RC and RCR model predictions against experimental pressure and flow data

using SimulationSetup

include("../theme.jl")

include(pkgdir(SimulationSetup, "..", "scripts", "validation", "pressure_model",
               "windkessel_model.jl"))

save_fig = true

set_theme!(my_thesis_theme)

fig1 = Figure(size=(800, 600) .* 0.5)
fig2 = Figure(size=(800, 600) .* 0.5)
ax_flow = Axis(fig1[1, 1], xlabel="Time [s]", ylabel="Flow [ml/s]",
               limits=((9T, 10T), nothing))
ax_pressure = Axis(fig2[1, 1], xlabel="Time [s]", ylabel="Pressure [mmHg]",
                   limits=((9T, 10T), nothing))

lines!(ax_pressure, sol_2_element.t, first.(sol_2_element.u), color=:blue, label="RC Model")
lines!(ax_pressure, sol_3_element.t, first.(sol_3_element.u), color=:red, label="RCR Model")
scatterlines!(ax_pressure, times_experiment .+ 9T, pressures_experiment; color=:black,
              marker=:circle, linestyle=:dash,
              markersize=10, label="Experiment")

axislegend(ax_pressure)

scatterlines!(ax_flow, times_experiment .+ 9T, flows_experiment, color=:black,
              linestyle=:dash, label="Experiment")
lines!(ax_flow, range(0, T, length=300) .+ 9T, t -> pulsatile_flow(t), color=:red,
       label="Fourier Series (fit)")
axislegend(ax_flow)

if save_fig
    dir = joinpath(fig_dir(), "windkessel_model")
    mkpath(dir)
    save(joinpath(dir, "flow_rates_experiment.pdf"), fig1)
    save(joinpath(dir, "pressures_experiment.pdf"), fig2)
else
    fig1
    fig2
end
