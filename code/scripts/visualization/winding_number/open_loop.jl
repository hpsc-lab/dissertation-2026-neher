using TrixiParticles
using CairoMakie

plot_open_curve = false

if plot_open_curve
    file = joinpath(data_dir(), "example_geometries", "inverted_open_curve.asc")
else
    file = joinpath(data_dir(), "example_geometries", "arbitrary_geometry.dxf")
end
geometry = load_geometry(file)

v = stack(geometry.vertices)
min_v = minimum(v, dims=2)
max_v = maximum(v, dims=2)

v = (v .- min_v) ./ (max_v .- min_v)
v[1, :] = v[1, :] .* 400 .+ 50
v[2, :] = v[2, :] .* 200 .+ 25

scaled_geometry = TrixiParticles.Polygon(v)

point_in_geometry_algorithm = WindingNumberJacobson(; geometry=scaled_geometry,
                                                    winding_number_factor=0.4,
                                                    hierarchical_winding=true)

grid = [SVector(x, y) for y in 1:250 for x in 1:500]

_, w = point_in_geometry_algorithm(scaled_geometry, grid; store_winding_number=true)

fig = Figure(size=(550, 275) .* 0.9)
ax = Axis(fig[1, 1])
ax.aspect = DataAspect()
hidedecorations!(ax)

hm = heatmap!(ax, 1:500, 1:250, reshape(w, 500, 250), colormap=:coolwarm,
              colorrange=plot_open_curve ? (-1, 1) : (-2, 2))

# Zeichne die Geometrie als Linie über die Heatmap
lines!(ax, v[1, :], v[2, :], color=:black, linewidth=3)

if !plot_open_curve
    let
        num_arrows = 11
        n = size(v, 2)
        step = max(1, n ÷ (num_arrows + 1))

        ps_x = Float64[]
        ps_y = Float64[]
        dirs_x = Float64[]
        dirs_y = Float64[]

        for i in 1:num_arrows
            idx = i * step
            if idx < n
                p1 = v[:, idx]
                p2 = v[:, idx + 1]
                push!(ps_x, p1[1])
                push!(ps_y, p1[2])
                push!(dirs_x, p2[1] - p1[1])
                push!(dirs_y, p2[2] - p1[2])
            end
        end
        arrows!(ax, ps_x, ps_y, dirs_x, dirs_y, color=:black, tipwidth=15, tiplength=8)
    end
end

cb = Colorbar(fig[1, 2], hm)
colsize!(fig.layout, 1, Aspect(1, 500 / 250))
resize_to_layout!(fig)

dir = joinpath(fig_dir(), "preprocessing", "sampling")
mkpath(dir)

if plot_open_curve
    save(joinpath(dir, "open_loop.png"), fig)
else
    save(joinpath(dir, "arbitrary_geometry.png"), fig)
end
fig
