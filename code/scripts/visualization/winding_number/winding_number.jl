using TrixiParticles
using CairoMakie

triangle = [125.0 375.0 250.0 125.0;
            175.0 175.0 350.0 175.0]

# Delete all edges but one
edge1 = deleteat!(TrixiParticles.Polygon(triangle), [2, 3])
edge2 = deleteat!(TrixiParticles.Polygon(triangle), [1, 3])
edge3 = deleteat!(TrixiParticles.Polygon(triangle), [1, 2])

algorithm = WindingNumberJacobson(; hierarchical_winding=false)

grid = [SVector(x, y) for y in 1:500 for x in 1:500]

_, w1 = algorithm(edge1, grid; store_winding_number=true)
_, w2 = algorithm(edge2, grid; store_winding_number=true)
_, w3 = algorithm(edge3, grid; store_winding_number=true)

w = w1 + w2 + w3

function save_heatmap(data, filename; clims=nothing)
    fig = Figure(size=(550, 500) .* 0.4)
    ax = Axis(fig[1, 1])
    ax.aspect = DataAspect()
    hidedecorations!(ax)

    if isnothing(clims)
        hm = heatmap!(ax, 1:500, 1:500, reshape(data, 500, 500), colormap=:coolwarm)
    else
        hm = heatmap!(ax, 1:500, 1:500, reshape(data, 500, 500), colormap=:coolwarm,
                      colorrange=clims)
    end
    cb = Colorbar(fig[1, 2], hm)
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)

    save(filename, fig)
end

dir = joinpath(fig_dir(), "preprocessing", "sampling")
mkpath(dir)

save_heatmap(w1, joinpath(dir, "triangle_w1.png"), clims=(-1, 1))
save_heatmap(w2, joinpath(dir, "triangle_w2.png"), clims=(-1, 1))
save_heatmap(w3, joinpath(dir, "triangle_w3.png"), clims=(-1, 1))
save_heatmap(w, joinpath(dir, "triangle_w.png"), clims=(-1, 1))
