# Plot pressure dynamics in healthy, stenosed, and arteriosclerotic vessels
# Includes validation/pressure_model/vessel_model.jl to generate sol, sol_stenosis, sol_arteriosclerosis
# Visualizes the effect of vessel pathologies (stenosis and arteriosclerosis) on pressure waveforms

using SimulationSetup
include("../theme.jl")

include(pkgdir(SimulationSetup, "..", "scripts", "validation", "pressure_model",
               "vessel_model.jl"))

save_fig = true

set_theme!(my_thesis_theme)

fig1 = Figure(size=(1600, 600) .* 0.5)
g = fig1[1, 1] = GridLayout()

# --- swapped axes: main = pressure (left), small = flow (right, outside main axis) ---
# left: main (large) pressure axis
ax_pressure = Axis(g[1, 1], xlabel="Time [s]", ylabel="Δp [mmHg]",
                   limits=((0, 7T), nothing))

# right: nested layout for legend (top) and small flow plot (bottom/outside main axis)
right = g[1, 2] = GridLayout()
legend_cell = right[1, 1]                       # top: legend
ax_inset_flow = Axis(right[2, 1], xlabel="Time [s]", ylabel="q [ml/s]")

# plot pressure traces on main axis (use your sol variables)
lines!(ax_pressure, sol.t, first.(sol.u), color=:blue, label="Physiological (healthy)")
lines!(ax_pressure, sol_stenosis.t, first.(sol_stenosis.u), color=:green,
       label="Stenotic vasoconstriction")
lines!(ax_pressure, sol_arteriosclerosis.t, first.(sol_arteriosclerosis.u), color=:red,
       label="Arteriosclerosis")
ylims!(ax_pressure, low=0)

# legend placed in the top cell of the right column (for the main pressure plot)
leg = Legend(legend_cell, ax_pressure)
leg.tellheight = true

# small flow trace in the right column below the legend (outside main axis)
t_flow = 0:dt:T
lines!(ax_inset_flow, t_flow, pulsatile_velocity_sin.(t_flow), color=:black)
xlims!(ax_inset_flow, 0, T)
ylims!(ax_inset_flow, low=0)

# optional: adjust spacing
colgap!(g, 10)
rowgap!(g, 8)

if save_fig
    dir = joinpath(fig_dir(), "windkessel_model")
    mkpath(dir)
    save(joinpath(dir, "pressures_vessel_model.pdf"), fig1)
else
    fig1
end
