using SimulationSetup
using CairoMakie
using Bessels
using CSV, DataFrames

save_fig = true

particle_spacing_factor = 30
pipe_radius = 0.0005

fluid_density = 1000.0
reynolds_number = 50
pressure_drop = 0.1
dynamic_viscosity = sqrt(fluid_density * pipe_radius^2 * pressure_drop /
                         (4 * reynolds_number))

# Analytical velocity evolution given in eq. 18 (Zhang et al., 2025)
function womersley_velocity_profile(r, t)
    omega = 1 # ω = 2π/T
    kinematic_viscosity = dynamic_viscosity / fluid_density
    alpha = pipe_radius * sqrt(omega / kinematic_viscosity)

    # pressure gradient magnitute for the frequency
    pressure_gradient = -25

    v_x = 0.0

    for n in [1]
        amp = im * pressure_gradient / (fluid_density * n * omega)

        term_1 = besselj0(alpha * sqrt(n) * im^(3 / 2) * r / pipe_radius)
        term_2 = besselj0(alpha * sqrt(n) * im^(3 / 2))

        exp_term = exp(im * n * omega * t)

        v_x += real(amp * (1 - term_1 / term_2) * exp_term)
    end

    return v_x
end

# Load results
output_directory = joinpath(out_dir(), "validation",
                            "open_boundaries", "pulsatile_channel_flow_3d")
data = CSV.read(joinpath(output_directory, "result_vx_dp_30.csv"), DataFrame)

times = round.(data[!, "time"], digits=2)
times_ref = round.(range(2pi, 4pi, step=0.51), digits=2)

r = range(-pipe_radius, pipe_radius, length=100)

data_indices = findall(t -> t in times_ref, times)
v_x_num = stack([eval(Meta.parse(str)) for str in data[!, "v_x_fluid_1"]][data_indices])  # 100 x Nt

# Analytisch: 100 x Nt
v_x_ana = stack([womersley_velocity_profile.(r, t) for t in times_ref])

Nt = length(times_ref)

# X-Offsets (Zeitachse “in Streifen”)
x_offsets = range(0, Nt * 5e-3; length=Nt)

# Alles um Offsets verschieben
for j in 1:Nt
    v_x_ana[:, j] .+= x_offsets[j]
    v_x_num[:, j] .+= x_offsets[j]
end

# Downsample numerische Punkte (wie 1:3:100)
idx_num = 1:3:100
r_num = r[idx_num]

# Ticks wie im Beispiel
xtick_pos = range(first(x_offsets), last(x_offsets), length=3)
xtick_lab = ["2π", "3π", "4π"]
ytick_pos = [-pipe_radius, 0.0, pipe_radius]
ytick_lab = ["R", "0", "R"]

# ---- Plot ----
fig = Figure(size=(1000, 300) .* 0.9)
ax = Axis(fig[1, 1];
          xlabel="time (s)",
          ylabel="radial coordinate",
          xticks=(collect(xtick_pos), xtick_lab),
          yticks=(ytick_pos, ytick_lab),)

# Vertikale Linien an jedem Zeit-Offset
vlines!(ax, collect(x_offsets); color=(:gray, 0.35), linewidth=1)

# Analytisch (gestrichelt) + numerisch (rot) pro Zeitscheibe
for j in 1:Nt
    scatter!(ax, v_x_num[idx_num, j], r_num; marker=:x, color=:red,# color=(:red, 0.5),
             markersize=8, label=(j == 1 ? "numerical" : nothing))
    lines!(ax, v_x_ana[:, j], r; color=:black, linestyle=:dash, linewidth=2,
           label=(j == 1 ? "analytical" : nothing))
end

xlims!(ax, first(x_offsets) - 1e-3, maximum(v_x_ana) + 1e-3)
ylims!(ax, -pipe_radius * 1.05, pipe_radius * 1.05)

Legend(fig[1, 2], ax)

if save_fig
    dir = joinpath(fig_dir(), "validation")
    mkpath(dir)
    save(joinpath(dir, "pulsatile_channel_flow_3d.pdf"), fig)
else
    display(fig)
end
