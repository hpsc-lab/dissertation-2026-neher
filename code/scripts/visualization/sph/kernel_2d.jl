include("../theme.jl")
using SimulationSetup
using TrixiParticles
using Random

save_fig = true

domain_size = 3.0
particle_spacing = 0.4
smoothing_length = 1.3 * particle_spacing
smoothing_kernel = SchoenbergCubicSplineKernel{2}()
compact_support = TrixiParticles.compact_support(smoothing_kernel, smoothing_length)

fig = Figure(size=(800, 600) .* 0.5)#my_figure_size(TEXT_WIDTH_PX; aspect_ratio=1.618))

ax = Axis3(fig[1, 1], elevation=0.06 * pi, # title="$(typeof(smoothing_kernel))", #aspect=:data,
           azimuth=0.2 * pi)

# ---- Particles ----
Random.seed!(1)

x_min, x_max = -domain_size / 2, domain_size / 2
y_min, y_max = -domain_size / 2, domain_size / 2
nx,
ny = round(Int, domain_size / particle_spacing),
     round(Int, domain_size / particle_spacing)
xg = LinRange(x_min, x_max, nx)
yg = LinRange(y_min, y_max, ny)

px = repeat(xg, inner=ny)
py = repeat(yg, outer=nx)
spacing = minimum(diff(xg))
perturb = 0.3 * spacing
px .= px .+ (rand(length(px)) .- 0.5) .* perturb
py .= py .+ (rand(length(py)) .- 0.5) .* perturb
pz = zeros(length(px))

rvals = [hypot(x, y) for (x, y) in zip(px, py)]
colors = [r <= compact_support ? :red : :gray for r in rvals]

scatter!(ax, px, py, pz; color=colors, markersize=8)

# ---- Kernel ----
xs = LinRange(-compact_support, compact_support, 100)
ys = LinRange(-compact_support, compact_support, 100)

Z = [TrixiParticles.kernel(smoothing_kernel, sqrt(x^2 + y^2),
                           smoothing_length)
     for x in xs, y in ys]
Z[Z .== 0.0] .= NaN

surface!(ax, xs, ys, Z, colormap=:coolwarm, transparency=true, shading=true, alpha=0.4,
         specular=0.75)

theta = LinRange(0, 2pi, 360)
xc = compact_support .* cos.(theta)
yc = compact_support .* sin.(theta)
zc = zeros(length(theta))
lines!(ax, xc, yc, zc, color=:black, linewidth=2, linestyle=:solid)

# --- center particle ----
scatter!(ax, [0.0], [0.0], [0.0]; color=:black, markersize=8)

# ---- vertical line ----
z_max = TrixiParticles.kernel(smoothing_kernel, 0, smoothing_length)
x0, y0 = 0.0, 0.0
z_top = z_max * 1.25 +
        0.75 * exp(-((x0 - 0.05)^2 + (y0 - 0.1)^2) / (2 * 0.6^2)) +
        0.2 * sin(pi * x0) * cos(pi * y0)
lines!(ax, [x0, x0], [y0, y0], [0.0, z_top]; color=:black, linewidth=2, linestyle=:dash)

# ---- field ----
xs = LinRange(-domain_size / 2, domain_size / 2, 200)
ys = LinRange(-domain_size / 2, domain_size / 2, 200)

# Basis: central Gaussian bell + light waves to add structure
Z = [z_max * 1.25 +
     0.75 * exp(-((x - 0.05)^2 + (y - 0.1)^2) / (2 * 0.6^2)) +   # central elevation
     0.2 * sin(pi * x) * cos(pi * y)                       # subtle ripples
     for x in xs, y in ys]

surface!(ax, xs, ys, Z; colormap=:viridis, shading=true,
         transparency=true, alpha=0.25, specular=0.75)

xs = LinRange(-domain_size / 2, domain_size / 2, 20)
ys = LinRange(-domain_size / 2, domain_size / 2, 20)
Z = [z_max * 1.25 +
     0.75 * exp(-((x - 0.05)^2 + (y - 0.1)^2) / (2 * 0.6^2)) +
     0.2 * sin(pi * x) * cos(pi * y)
     for x in xs, y in ys]
wireframe!(ax, xs, ys, Z; color=:gray, linewidth=0.5)

scatter!(ax, [x0], [y0], [z_top]; color=:black, markersize=12, marker=:xcross)
zlims!(ax, low=0)

ax.xlabel = "x"
ax.ylabel = "y"
ax.zlabel = "f"
# ax.zticks = ([], [])
text!(fig[1, 1], "W(r-r', h)"; position=(0.0, 0.4, 0.65z_max), #rotation=pi / 2,
      #   align=(:left, :left),
      fontsize=14, color=:black)
text!(fig[1, 1], "f(r)"; position=(-1.0, 1.2, 1.3z_max), #rotation=pi / 2,
      #   align=(:left, :left),
      fontsize=14, color=:black)

if save_fig
    dir = joinpath(fig_dir(), "sph")
    mkpath(dir)
    save(joinpath(dir, "kernel_2d.png"), fig)
else
    fig
end
