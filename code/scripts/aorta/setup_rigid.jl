using SimulationSetup
using TrixiParticles, OrdinaryDiffEqLowStorageRK
# using ThreadPinning
# pinthreads(:numa)

# ==========================================================================================
# ==== Configuration

# Parse command-line arguments
param1, param2,
param3 = length(ARGS) >= 3 ? (ARGS[1], parse(Float64, ARGS[2]), ARGS[3]) :
         ("F10", 0.001, "1.0.21")
scenario_override = length(ARGS) >= 4 ? Symbol(ARGS[4]) : :normotensive

# Configure simulation with command-line parameters
version = VersionNumber(param3)
scenario = isnothing(scenario_override) ? infer_scenario(version) : scenario_override

set_config!(version=version, scenario=scenario)
initialize_code_version!()

# Patient identifier
subject = param1
# Resolution
particle_spacing = param2

wall_thickness = 4 * particle_spacing

trixi_include(@__MODULE__, joinpath(@__DIR__, "setup_fluid.jl"), subject=subject,
              factor_buffer_size=2, particle_spacing=particle_spacing, coord_eltype=Float64)

print_startup_message(param_sim; simulation_setup=:rigid,
                      flow_rate_correction_factor=flow_rate_correction_factor,
                      factor_buffer_size=factor_buffer_size, sound_speed=sound_speed)

# Prepare the vessel wall particles.
# The pressure field in "aorta_boundary" is repurposed to store the signed distances.
signed_distances_wall = ic_dict["aorta_boundary"].pressure

# Select particles within the vessel wall thickness
keep_indices = signed_distances_wall .<= wall_thickness

# Extract coordinates of wall particles and initialize the vessel wall
boundary_coordinates = ic_dict["aorta_boundary"].coordinates[:, keep_indices]
vessel_wall = InitialCondition(; coordinates=boundary_coordinates, density=density_blood,
                               particle_spacing)

# ==========================================================================================
# ==== Boundary
boundary_model = BoundaryModelDummyParticles(vessel_wall.density, vessel_wall.mass,
                                             state_equation=state_equation,
                                             AdamiPressureExtrapolation(),
                                             smoothing_kernel, smoothing_length,
                                             viscosity=viscosity)

boundary_system = WallBoundarySystem(vessel_wall, boundary_model)

# ==========================================================================================
# ==== Simulation
tspan = (0.0, ncycles() * T)

min_corner = minimum(vessel_wall.coordinates .- 2 * particle_spacing, dims=2)
max_corner = maximum(vessel_wall.coordinates .+ 2 * particle_spacing, dims=2)

nhs = GridNeighborhoodSearch{3}(; cell_list=FullGridCellList(; min_corner, max_corner),
                                update_strategy=ParallelUpdate())

semi = Semidiscretization(fluid_system, open_boundary,
                          boundary_system, neighborhood_search=nhs,
                          parallelization_backend=PolyesterBackend())

ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)

output_directory = joinpath(out_dir(), "out_$(current_scenario())", "aorta",
                            "$subject", "rigid", "dp_$(particle_spacing)")

pp_callback = PostprocessCallback(; dt=0.01, output_directory,
                                  filename="resulting_pressures", write_csv=true,
                                  write_file_interval=1, pp_functions(boundaries)...)

saving_callback = SolutionSavingCallback(dt=0.01, prefix="", overwrite=true,
                                         output_directory=output_directory)

extra_callback = SortingCallback(; interval=1000)

callbacks = CallbackSet(info_callback, extra_callback, saving_callback, pp_callback,
                        UpdateCallback())

sol = solve(ode, RDPK3SpFSAL35(),
            abstol=1e-7, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=1e-4, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            dtmax=1e-2, # Limit stepsize to prevent crashing
            save_everystep=false, callback=callbacks);
