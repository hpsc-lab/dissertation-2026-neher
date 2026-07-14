# include("../theme.jl")
using SimulationSetup
using TrixiParticles
using CairoMakie

save_fig = true

n_dims = 2

kernel_func(kernel, r, h) = TrixiParticles.kernel(kernel, r, h)
kernel_deriv_func(kernel, r, h) = TrixiParticles.kernel_deriv(kernel, r, h)

kernel = WendlandC2Kernel{n_dims}()

# Parameters
h = 1.0
r_s = TrixiParticles.compact_support(kernel, h)
r = range(0, stop=r_s, length=600)

kernel_max = kernel_func(kernel, 0, h)

y = [kernel_func(kernel, ri, h) / kernel_max for ri in r]

fig = Figure(size=(800, 400) .* 0.5)
ax = Axis(fig[1, 1], xlabel="ζ", ylabel="w_δv")
lines!(ax, r, y, color=:blue, label="smooth activation")
lines!(ax, [0, 0], [0, 1], color=:red)  # rote vertikale Linie bei x=0 von y=0 bis y=1
lines!(ax, [-0.5, 0.0], [1.0, 1.0], color=:red, label="discontinuous activation")  # rote vertikale Linie bei x=0 von y=0 bis y=1
lines!(ax, [0, r_s * 2], [0, 0], color=:red)  # rote horizontale Linie bei y=0 ab x=0
axislegend(ax, position=:rt)
ylims!(ax, low=-0.1, high=1.1)
xlims!(ax, low=-0.1, high=2.1)

if save_fig
    dir = joinpath(fig_dir(), "open_boundaries")
    mkpath(dir)
    save(joinpath(dir, "kernel_weight.pdf"), fig)
else
    fig
end
