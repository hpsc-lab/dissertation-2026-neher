using SimulationSetup
using TrixiParticles, OrdinaryDiffEqLowStorageRK
# using ThreadPinning
# pinthreads(:numa)

# ==========================================================================================
# ==== Configuration

# Use command-line arguments if provided; otherwise fall back to the default values
param1, param2,
param3 = length(ARGS) >= 3 ? (ARGS[1], parse(Float64, ARGS[2]), ARGS[3]) :
         ("F10", 0.0005, "1.0.21")
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

trixi_include(@__MODULE__, joinpath(@__DIR__, "setup_fluid.jl"), subject=subject,
              factor_buffer_size=4, particle_spacing=particle_spacing, coord_eltype=Float64)

wall_thickness = 0.002
vessel_thickness = 0.002
density_vessel = 1000.0

E_modul_vessel = param_sim.subject_parameters.E
E_modul_tissue = E_modul_vessel
nu_vessel = param_sim.subject_parameters.nu
nu_tissue = 0.0 # TODO: explain this.

print_startup_message(param_sim; simulation_setup=:elastic, wall_thickness=wall_thickness,
                      vessel_thickness=vessel_thickness,
                      flow_rate_correction_factor=flow_rate_correction_factor,
                      factor_buffer_size=factor_buffer_size, sound_speed=sound_speed)

E_moduli = fill(E_modul_tissue, nparticles(ic_dict["aorta_boundary"]))
poisson_ratios = fill(nu_tissue, nparticles(ic_dict["aorta_boundary"]))

# Prepare the vessel wall particles.
# The pressure field in "aorta_boundary" is repurposed to store the signed distances.
signed_distances_wall = ic_dict["aorta_boundary"].pressure

# Select particles within the vessel wall thickness
keep_indices = signed_distances_wall .<= wall_thickness
vessel_indices = signed_distances_wall .<= vessel_thickness

E_moduli[vessel_indices] .= E_modul_vessel
poisson_ratios[vessel_indices] .= nu_vessel

# Extract coordinates of wall particles and initialize the vessel wall
boundary_coordinates = ic_dict["aorta_boundary"].coordinates[:, keep_indices]
vessel_wall = InitialCondition(; coordinates=boundary_coordinates, density=density_vessel,
                               particle_spacing)
deleteat!(E_moduli, .!keep_indices)
deleteat!(poisson_ratios, .!keep_indices)

# ==========================================================================================
# ==== Structure
clamped_candidates = Int[]
for key in keys(boundaries)
    bz = boundary_zones[key]
    geometry_planar = boundaries[key].geometry
    geometry = extrude_geometry(geometry_planar, bz.zone_width)

    signed_distance_field = SignedDistanceField(geometry, particle_spacing;
                                                use_for_boundary_packing=true,
                                                max_signed_distance=wall_thickness)

    boundary_sampled = sample_boundary(signed_distance_field;
                                       boundary_density=density_vessel,
                                       boundary_thickness=wall_thickness)

    coords_rel = boundary_sampled.coordinates .- bz.zone_origin .-
                 3 * particle_spacing .* bz.face_normal
    coords_rel = reinterpret(reshape, SVector{3, eltype(coords_rel)}, coords_rel)

    keep = findall(x -> signbit(TrixiParticles.dot(x, bz.face_normal)), coords_rel)

    ids = TrixiParticles.find_too_close_particles(vessel_wall.coordinates,
                                                  boundary_sampled.coordinates[:, keep],
                                                  2 * particle_spacing)

    append!(clamped_candidates, ids)
end

hydrodynamic_densites = density_blood * ones(size(vessel_wall.density))
hydrodynamic_masses = hydrodynamic_densites * particle_spacing^ndims(vessel_wall)

boundary_model = BoundaryModelDummyParticles(hydrodynamic_densites, hydrodynamic_masses,
                                             AdamiPressureExtrapolation(),
                                             smoothing_kernel, smoothing_length,
                                             viscosity=viscosity)

structure_system = TotalLagrangianSPHSystem(vessel_wall; smoothing_kernel,
                                            smoothing_length, young_modulus=E_moduli,
                                            poisson_ratio=poisson_ratios,
                                            boundary_model=boundary_model,
                                            penalty_force=PenaltyForceGanzenmueller(alpha=0.05),
                                            viscosity=ArtificialViscosityMonaghan(alpha=0.05),
                                            clamped_particles=unique!(clamped_candidates))

# ==========================================================================================
# ==== Simulation
tspan = (0.0, ncycles() * T)

min_corner = minimum(vessel_wall.coordinates .- 10 * vessel_thickness, dims=2)
max_corner = maximum(vessel_wall.coordinates .+ 10 * vessel_thickness, dims=2)

nhs = GridNeighborhoodSearch{3}(; cell_list=FullGridCellList(; min_corner, max_corner),
                                update_strategy=ParallelUpdate())

semi = Semidiscretization(fluid_system, open_boundary, structure_system,
                          neighborhood_search=nhs,
                          parallelization_backend=PolyesterBackend())

ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)

output_directory = joinpath(out_dir(), "out_$(current_scenario())", "aorta", "$subject",
                            "elastic", "dp_$(particle_spacing)_t_$(vessel_thickness)")

pp_callback = PostprocessCallback(; dt=0.01, output_directory,
                                  filename="resulting_pressures", write_csv=true,
                                  write_file_interval=1, pp_functions(boundaries)...)

saving_callback = SolutionSavingCallback(dt=0.01, prefix="", overwrite=true,
                                         output_directory=output_directory)

extra_callback = SortingCallback(; interval=1000)

callbacks = CallbackSet(info_callback, extra_callback, saving_callback, pp_callback,
                        UpdateCallback())

sol = solve(ode, RDPK3SpFSAL35(),
            maxiters=300_000,
            abstol=1e-7, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=1e-4, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            dtmax=1e-2, # Limit stepsize to prevent crashing
            save_everystep=false, callback=callbacks);
