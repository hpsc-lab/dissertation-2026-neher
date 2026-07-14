# Wall Shear Stress (WSS) Validation for 2D Poiseuille Flow

# This script validates the turbulence model's wall shear stress calculation against
# analytical solutions for a 2D Poiseuille channel flow (constant pressure gradient).

# The plot shows:
# - Relative error in WSS calculation vs. analytical solution
# - Two computation strategies compared:
#   1. WSS on boundary particles (red line)
#   2. WSS extrapolated to physical wall (blue line)
# - X-axis: Particle resolution (particles per wall distance)
# - Y-axis: Relative error [%]
using SimulationSetup
using TrixiParticles
using CairoMakie

include("../theme.jl")
set_theme!(my_thesis_theme)

save_fig = true

on_wall = parse(Bool, ARGS[1])

# Analytical velocity evolution given in eq. 16 (Zhang et al., 2025)
function poiseuille_velocity(y, t)

    # Base profile (stationary part)
    base_profile = (pressure_drop / (2 * dynamic_viscosity * flow_length)) * y *
                   (y - wall_distance)

    # Transient terms (Fourier series)
    transient_sum = 0.0

    for n in 0:10  # Limit to 10 terms for convergence
        coefficient = (4 * pressure_drop * wall_distance^2) /
                      (dynamic_viscosity * flow_length * pi^3 * (2 * n + 1)^3)

        sine_term = sin(pi * y * (2 * n + 1) / wall_distance)

        exp_term = exp(-((2 * n + 1)^2 * pi^2 * dynamic_viscosity * t) /
                       (fluid_density * wall_distance^2))

        transient_sum += coefficient * sine_term * exp_term
    end

    # Total velocity
    v_x = base_profile + transient_sum

    return -v_x
end

function compute_wss(; divisors, on_wall=true)
    wall_distance = 0.001 # distance between top and bottom wall
    flow_length = 0.004 # distance between inflow and outflow

    # Make sure that the kernel support of fluid particles at a boundary is always fully sampled
    boundary_layers = 4

    domain_size = (flow_length, wall_distance)

    fluid_density = 1000.0
    reynolds_number = 50
    pressure_drop = 0.1
    dynamic_viscosity = sqrt(fluid_density * wall_distance^3 * pressure_drop /
                             (8 * flow_length * reynolds_number))

    v_max = wall_distance^2 * pressure_drop / (8 * dynamic_viscosity * flow_length)

    WSS_calculated_physical_wall = Float64[]
    WSS_calculated = Float64[]

    for divisor in divisors
        @info divisor
        particle_spacing = wall_distance / divisor

        pipe = RectangularTank(particle_spacing, domain_size, domain_size, fluid_density,
                               velocity=(pos) -> (poiseuille_velocity(pos[2], 10.0), 0.0),
                               n_layers=boundary_layers, faces=(false, false, true, true))

        kinematic_viscosity = dynamic_viscosity / fluid_density
        viscosity = ViscosityAdami(nu=kinematic_viscosity)

        state_equation = StateEquationCole(; sound_speed=10,
                                           reference_density=fluid_density,
                                           exponent=1)

        smoothing_kernel = SchoenbergQuinticSplineKernel{2}()
        smoothing_length = 1.1 * particle_spacing
        fluid_system = WeaklyCompressibleSPHSystem(pipe.fluid;
                                                   density_calculator=ContinuityDensity(),
                                                   state_equation, smoothing_kernel,
                                                   smoothing_length, viscosity=viscosity)
        wall = pipe.boundary
        boundary_model = BoundaryModelDummyParticles(wall.density, wall.mass,
                                                     AdamiPressureExtrapolation(),
                                                     state_equation=state_equation,
                                                     viscosity=viscosity,
                                                     smoothing_kernel, smoothing_length)

        boundary_system = WallBoundarySystem(wall, boundary_model)

        semi = Semidiscretization(fluid_system, boundary_system)

        ode = semidiscretize(semi, (0.0, 1.0))

        v_ode, u_ode = ode.u0.x
        on_wall && TrixiParticles.update_systems_and_nhs(v_ode, u_ode, semi, 0.0)

        isotropic_constant = 0.0
        # isotropic_constant = 6.6e-3
        smagorinsky_constant = 0.0
        # smagorinsky_constant = 0.12
        dynamic_viscosity_ = dynamic_viscosity
        # dynamic_viscosity_ = 0.0

        turbulence_model_fluid = SPSTurbulenceModelDalrymple(fluid_system.initial_condition;
                                                             smallest_length_scale=particle_spacing,
                                                             dynamic_viscosity=dynamic_viscosity_,
                                                             smagorinsky_constant,
                                                             isotropic_constant)

        n_points = round(Int, (flow_length / 50) / particle_spacing)
        sample_points_physical_wall = RectangularShape(particle_spacing, (n_points, 1),
                                                       (flow_length / 2, wall_distance),
                                                       density=fluid_density,
                                                       place_on_shell=true)

        turbulence_model_wall = SPSTurbulenceModelDalrymple(sample_points_physical_wall;
                                                            smallest_length_scale=particle_spacing,
                                                            dynamic_viscosity=dynamic_viscosity_,
                                                            smagorinsky_constant,
                                                            isotropic_constant,
                                                            only_wall_shear_stress=true)

        sample_points_boundary = RectangularShape(particle_spacing, (n_points, 1),
                                                  (flow_length / 2, wall_distance),
                                                  density=fluid_density,
                                                  place_on_shell=false)

        turbulence_model_boundary = SPSTurbulenceModelDalrymple(sample_points_boundary;
                                                                smallest_length_scale=particle_spacing,
                                                                dynamic_viscosity=dynamic_viscosity,
                                                                smagorinsky_constant=0.0,
                                                                isotropic_constant=0.0,
                                                                only_wall_shear_stress=true)

        SimulationSetup.calculate_fluid_stress_tensor!(fluid_system, turbulence_model_fluid,
                                                       v_ode, u_ode, semi)

        SimulationSetup.calculate_wall_shear_stress!(turbulence_model_wall,
                                                     turbulence_model_fluid, fluid_system,
                                                     v_ode, u_ode, semi)

        SimulationSetup.calculate_wall_shear_stress!(turbulence_model_boundary,
                                                     turbulence_model_fluid, fluid_system,
                                                     v_ode, u_ode, semi)

        push!(WSS_calculated_physical_wall,
              sum(TrixiParticles.norm.(turbulence_model_wall.field_variables.stress_vectors)) /
              n_points)

        push!(WSS_calculated,
              sum(TrixiParticles.norm.(turbulence_model_boundary.field_variables.stress_vectors)) /
              n_points)
    end

    return WSS_calculated, WSS_calculated_physical_wall
end

wall_distance = 0.001 # distance between top and bottom wall
flow_length = 0.004 # distance between inflow and outflow

# Make sure that the kernel support of fluid particles at a boundary is always fully sampled
boundary_layers = 4

domain_size = (flow_length, wall_distance)

fluid_density = 1000.0
reynolds_number = 50
pressure_drop = 0.1
dynamic_viscosity = sqrt(fluid_density * wall_distance^3 * pressure_drop /
                         (8 * flow_length * reynolds_number))

v_max = wall_distance^2 * pressure_drop / (8 * dynamic_viscosity * flow_length)

x_values = [25, 50, 100, 200, 400, 800]
WSS_calculated, WSS_calculated_physical_wall = compute_wss(; divisors=x_values, on_wall)

WSS_analytic = 4 * dynamic_viscosity * v_max / wall_distance * ones(length(WSS_calculated))

error_WSS = abs.(WSS_calculated .- WSS_analytic) ./ abs.(WSS_analytic) .* 100
error_WSS_physical_wall = abs.(WSS_calculated_physical_wall .- WSS_analytic) ./
                          abs.(WSS_analytic) .* 100

fig = Figure(size=(800, 600) .* 0.5)
ax = Axis(fig[1, 1],
          xlabel="particles per wall distance",
          ylabel="rel. error [%]")

scatterlines!(ax, x_values, error_WSS, markersize=10, label="on boundary particle",
              color=:red)
scatterlines!(ax, x_values, error_WSS_physical_wall, markersize=10,
              label="on physical wall", color=:blue)
xlims!(ax, low=0)
ylims!(ax, low=0)
if on_wall
    axislegend(ax, position=:rt)
else
    axislegend(ax, position=:rb)
end

if save_fig
    save(joinpath(fig_dir(), "validation",
                  "wall_shear_stress$(on_wall ? "_on_wall" : "").pdf"), fig)
else
    fig
end
