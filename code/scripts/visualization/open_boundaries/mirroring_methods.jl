using CairoMakie
using SimulationSetup
using TrixiParticles

save_fig = true

# Create a fluid domain with open boundaries and apply a mirroring method
function mirror(pressure_function, mirror_method;
                particle_spacing=0.05, domain_size=(2.0, 1.0))
    # Setup fluid domain with rectangular shape
    domain_fluid = RectangularShape(particle_spacing,
                                    round.(Int, domain_size ./ particle_spacing),
                                    (0.0, 0.0), density=1000.0,
                                    pressure=pressure_function)

    smoothing_length = 1.2 * particle_spacing
    smoothing_kernel = WendlandC2Kernel{2}()
    fluid_system = EntropicallyDampedSPHSystem(domain_fluid; smoothing_kernel,
                                               smoothing_length, sound_speed=1.0)

    fluid_system.cache.density .= domain_fluid.density

    # Setup outflow boundary at x = domain_size[1]
    plane_out = ([domain_size[1], 0.0], [domain_size[1], domain_size[2]])

    outflow = BoundaryZone(; boundary_face=plane_out, boundary_type=OutFlow(),
                           face_normal=[-1.0, 0.0],
                           open_boundary_layers=10, density=1000.0, particle_spacing)
    open_boundary_out = OpenBoundarySystem(outflow; fluid_system,
                                           boundary_model=BoundaryModelMirroringTafuni(),
                                           buffer_size=0)
    open_boundary_out.boundary_zone_indices .= 1

    semi = Semidiscretization(fluid_system, open_boundary_out)

    TrixiParticles.initialize_neighborhood_searches!(semi)

    v_open_boundary = zero(outflow.initial_condition.velocity)
    v_fluid = vcat(domain_fluid.velocity, domain_fluid.pressure')

    TrixiParticles.set_zero!(open_boundary_out.cache.pressure)

    TrixiParticles.extrapolate_values!(open_boundary_out, mirror_method,
                                       v_open_boundary, v_fluid,
                                       outflow.initial_condition.coordinates,
                                       domain_fluid.coordinates, semi)

    # Setup inflow boundary at x = 0.0
    plane_in = ([0.0, 0.0], [0.0, domain_size[2]])

    inflow = BoundaryZone(; boundary_face=plane_in, boundary_type=InFlow(),
                          face_normal=[1.0, 0.0],
                          open_boundary_layers=10, density=1000.0, particle_spacing)
    open_boundary_in = OpenBoundarySystem(inflow; fluid_system,
                                          boundary_model=BoundaryModelMirroringTafuni(),
                                          buffer_size=0)
    open_boundary_in.boundary_zone_indices .= 1

    semi = Semidiscretization(fluid_system, open_boundary_in)
    TrixiParticles.initialize_neighborhood_searches!(semi)

    v_open_boundary = zero(inflow.initial_condition.velocity)

    TrixiParticles.set_zero!(open_boundary_in.cache.pressure)

    TrixiParticles.extrapolate_values!(open_boundary_in, mirror_method,
                                       v_open_boundary, v_fluid,
                                       inflow.initial_condition.coordinates,
                                       domain_fluid.coordinates, semi)

    return fluid_system, open_boundary_in, open_boundary_out, v_fluid
end

# Interpolate pressure along a horizontal line through the domain
function interpolate_pressure(mirror_method, pressure_func; particle_spacing=0.05)
    fluid_system, open_boundary_in, open_boundary_out,
    v_fluid = mirror(pressure_func, mirror_method)

    p_fluid = [TrixiParticles.current_pressure(v_fluid, fluid_system, particle)
               for particle in TrixiParticles.each_active_particle(fluid_system)]

    fluid_system.initial_condition.pressure .= p_fluid
    open_boundary_in.initial_condition.pressure .= open_boundary_in.cache.pressure
    open_boundary_out.initial_condition.pressure .= open_boundary_out.cache.pressure

    entire_domain = union(fluid_system.initial_condition,
                          open_boundary_in.initial_condition,
                          open_boundary_out.initial_condition)

    smoothing_length = 1.2 * particle_spacing
    smoothing_kernel = WendlandC2Kernel{2}()

    # Use a fluid system to interpolate the pressure across entire domain
    interpolation_system = WeaklyCompressibleSPHSystem(entire_domain;
                                                       density_calculator=ContinuityDensity(),
                                                       state_equation=nothing,
                                                       smoothing_kernel, smoothing_length)
    interpolation_system.pressure .= entire_domain.pressure

    semi = Semidiscretization(interpolation_system)
    ode = semidiscretize(semi, (0, 0))
    v_ode, u_ode = ode.u0.x

    result = interpolate_line([-0.5, 0.5], [2.5, 0.5], 100, semi,
                              interpolation_system, v_ode, u_ode)

    return result.pressure
end

# Define test pressure functions
pressure_func_1(pos) = 2pos[1]
pressure_func_2(pos) = pos[1]^2
pressure_func_3(pos) = 1.2cos(2pi * pos[1]) + 0.5

pressure_funcs = [
    (pressure_func_1, "p(x) = 2x"),
    (pressure_func_2, "p(x) = x²"),
    (pressure_func_3, "p(x) = 1.2cos(2πx) + 0.5")
]

# Create figure
fig = Figure(resolution=(1200, 350) .* 0.8)

label_ = ["Simple", "1st Order", "0th Order"]
linestyles = [:solid, :dot, :dash]
color = [:red, :blue, :black]

# Store legend elements
legend_elements = []

# Calculate all pressure values first to find global limits
all_pressures = []
for (pressure_func, title) in pressure_funcs
    pressures = interpolate_pressure.([
                                          SimpleMirroring(),
                                          FirstOrderMirroring(),
                                          ZerothOrderMirroring()
                                      ],
                                      pressure_func)
    push!(all_pressures, pressures)
end

# Set global y-limits for consistent scaling across all plots
global_y_min = -1.2
global_y_max = 5.5

for (i, (pressure_func, title)) in enumerate(pressure_funcs)
    pressures = all_pressures[i]

    # Create subplot (only left plot has y-label)
    if i == 1
        ax = Axis(fig[1, i], xlabel="x", ylabel="p(x)", title=title)
    else
        ax = Axis(fig[1, i], xlabel="x", title=title)
        hideydecorations!(ax, label=false, ticklabels=true, ticks=true, grid=false)
    end

    # Set same y-limits for all axes
    ylims!(ax, global_y_min, global_y_max)

    # Gray background box for fluid domain (x ∈ [0, 2])
    x_min, x_max = 0.0, 2.0

    poly!([Point2f(x_min, global_y_min), Point2f(x_max, global_y_min),
              Point2f(x_max, global_y_max), Point2f(x_min, global_y_max)],
          color=(:lightgray, 0.25), strokewidth=0)

    # Plot pressure profiles
    x = collect(range(-0.5, stop=2.5, length=length(first(pressures))))

    for (j, p) in enumerate(pressures)
        l = lines!(ax, x, p, linewidth=2,
                   color=color[j], linestyle=linestyles[j])
        # Save line elements only from first plot for legend
        if i == 1
            push!(legend_elements, l)
        end
    end
end

# Move axes closer together
colgap!(fig.layout, 0)

# Add legend outside on the right
Legend(fig[1, 4], legend_elements, label_, framevisible=true)

if save_fig
    dir = joinpath(fig_dir(), "open_boundaries")
    mkpath(dir)
    save(joinpath(dir, "extrapolated_values_with_different_mirroring_methods.pdf"), fig)
else
    fig
end
