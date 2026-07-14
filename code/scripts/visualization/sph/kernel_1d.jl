include("../theme.jl")
using SimulationSetup
using TrixiParticles

save_fig = true

n_dims = 2

kernel_func(kernel, r, h) = TrixiParticles.kernel(kernel, r, h)
kernel_deriv_func(kernel, r, h) = TrixiParticles.kernel_deriv(kernel, r, h)

kernel_1 = SchoenbergCubicSplineKernel{n_dims}()
kernel_2 = SchoenbergQuinticSplineKernel{n_dims}()
kernel_3 = WendlandC2Kernel{n_dims}()
kernel_gauss = GaussianKernel{n_dims}()

# Parameters
h = 1.0
r = range(0, stop=3.8h, length=600)
x = r ./ h

kernels = [kernel_1, kernel_2, kernel_3, kernel_gauss]
titles = ["Cubic Spline", "Quintic Spline", "Wendland C2", "Gaussian"]

labels = ["W(q)", "dW(q)"]

fig = Figure(size=(1000, 400) .* 0.65)
gl = fig[1, 1] = GridLayout(tellwidth=false, tellheight=true)
axes = Axis[]

for (i, ttl) in enumerate(titles)
    ax = Axis(gl[1, i];
              title=ttl,
              xlabel="r/h",
              ylabel=i == 1 ? "f(r, h)" : "",
              yminorticksvisible=true,
              yminorticks=IntervalsBetween(5),
              yminorgridvisible=true,
              xminorticksvisible=true,
              spinewidth=1.0,
              xgridvisible=true,
              ygridvisible=true,
              xticksize=8,
              yticksize=8,
              yticksmirrored=false)

    # hlines!(ax, 0; linewidth=1, color=:black)

    push!(axes, ax)
end

# Gemeinsame Achsen koppeln
linkxaxes!(axes...)
linkyaxes!(axes...)

y_gauss = [kernel_func(kernel_gauss, ri, h) for ri in r]
dy_gauss = [kernel_deriv_func(kernel_gauss, ri, h) for ri in r]

for (i, kernel) in enumerate(kernels)
    y = [kernel_func(kernel, ri, h) for ri in r]
    dy = [kernel_deriv_func(kernel, ri, h) for ri in r]

    if i < length(kernels)
        lines!(axes[i], x, [iszero(v) ? NaN : v for v in y], color=:black, linewidth=3)
        lines!(axes[i], x, [iszero(v) ? NaN : v for v in dy], color=:red, linewidth=3)
        lines!(axes[i], x, y_gauss, color=:black, linestyle=:dot)
        lines!(axes[i], x, dy_gauss, color=:red, linestyle=:dot)
    else
        lines!(axes[i], x, y_gauss, color=:black, linewidth=3)
        lines!(axes[i], x, dy_gauss, color=:red, linewidth=3)
    end

    xlims!(axes[i], 0, last(x))
    ylims!(axes[i], -1, 1)

    i > 1 && hideydecorations!(axes[i], grid=false, minorgrid=false, minorticks=false)
end

colgap!(gl, 0)

if save_fig
    export_fig(joinpath(fig_dir(), "sph", "kernel_1d"), fig; save_pdf=true)
else
    fig
end
