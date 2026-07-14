# Preprocessing script for creating simulation-ready initial conditions.

# Converts preprocessed aorta geometries into simulation-ready particle distributions
# and boundary conditions using signed distance fields and particle packing algorithms.

# Steps:
#     1. Load aorta geometry and boundary surface files
#     2. Extrude boundary surfaces for transition zones
#     3. Merge and union all geometry components
#     4. Sample interior domain with particles using winding number algorithm
#     5. Sample boundary surfaces and apply boundary conditions
#     6. Output initial conditions as VTK files

# Command-line arguments:
#     param1 (default: "F10"): Subject identifier (F01-F16)
#     param2 (default: 0.001): Particle spacing in meters

using SimulationSetup
using TrixiParticles, OrdinaryDiffEqLowStorageRK
# using ThreadPinning
# pinthreads(:numa)

# Use command-line arguments if provided; otherwise fall back to the default values
param1, param2,
param3 = length(ARGS) >= 3 ? (ARGS[1], parse(Float64, ARGS[2]), ARGS[3]) :
         ("F10", 0.001, "1.0.1")

# Configure simulation version and auto-detect code version from Git
set_config!(version=VersionNumber(param3))
initialize_code_version!()

# ==========================================================================================
# ==== Load geometries
subject = param1
particle_spacing = param2
boundary_thickness = particle_spacing < 0.005 ? 6 * particle_spacing : 10 * particle_spacing
transition_length = 5e-3
boundary_zone_width = 14 * particle_spacing
extrusion_length = boundary_zone_width + transition_length

geometry_boundaries = boundary_names()

files = joinpath.(data_dir(), "aorta_preprocessed", "v$(current_version().major)",
                  subject, subject * "_" .* geometry_boundaries .* ".stl")
file_aorta = joinpath(data_dir(), "aorta_preprocessed", "v$(current_version().major)",
                      subject, subject * ".stl")

SimulationSetup.check_outlet_configuration!(files, geometry_boundaries)

output_directory = joinpath(data_dir(), "aorta_initial_condition",
                            "v$(current_version().major)", "packed_results_" * subject,
                            "dp_$(particle_spacing)")

# ==========================================================================================
# ==== Preprocess geometries
geometries = load_geometry.(pushfirst!(files, file_aorta))
geometries_extruded = extrude_geometry.(geometries[2:end], extrusion_length,
                                        omit_bottom_face=true)

geometry = union(first(geometries), geometries_extruded...)

# ==========================================================================================
# ==== Packing parameters
place_on_shell = false
pack_boundary = true

maxiters_ = 5000
abstol = 1e-7
reltol = 1e-4

h_factor = 1.0
h_factor_interpolation = 1.0
smoothing_kernel = SchoenbergQuinticSplineKernel{3}()

background_pressure = 1.0
boundary_compress_factor = 0.9
smoothing_length_interpolation = h_factor_interpolation * particle_spacing

# ==========================================================================================
# ==== Sample geometry
density = 1.0

signed_distance_field = SignedDistanceField(geometry, particle_spacing;
                                            use_for_boundary_packing=true,
                                            max_signed_distance=boundary_thickness)

point_in_geometry_algorithm = WindingNumberJacobson(; geometry)

# Returns `InitialCondition`
shape_sampled = ComplexShape(geometry; particle_spacing, density, grid_offset=0.0,
                             max_nparticles=10^10, point_in_geometry_algorithm)

shape_sampled.mass .= density * TrixiParticles.volume(geometry) /
                      nparticles(shape_sampled)

trixi2vtk(shape_sampled, output_directory=output_directory,
          filename="sampled_aorta")

boundary_sampled = sample_boundary(signed_distance_field; boundary_density=density,
                                   boundary_thickness, place_on_shell)

trixi2vtk(boundary_sampled, output_directory=output_directory,
          filename="sampled_aorta_boundary")

# ==========================================================================================
# ==== Packing
smoothing_length = h_factor * particle_spacing

packing_system = ParticlePackingSystem(shape_sampled;
                                       smoothing_kernel=smoothing_kernel,
                                       smoothing_length=smoothing_length,
                                       smoothing_length_interpolation=smoothing_length,
                                       signed_distance_field, place_on_shell=place_on_shell,
                                       background_pressure)
boundary_system = ParticlePackingSystem(boundary_sampled;
                                        smoothing_kernel=smoothing_kernel,
                                        smoothing_length=smoothing_length,
                                        smoothing_length_interpolation=smoothing_length,
                                        boundary_compress_factor=boundary_compress_factor,
                                        is_boundary=true, signed_distance_field,
                                        place_on_shell=place_on_shell, background_pressure)

# ==========================================================================================
# ==== Simulation
semi = Semidiscretization(packing_system, boundary_system)

# Use a high `tspan` to guarantee that the simulation runs at least for `maxiters`
tspan = (0, 10000.0)
ode = semidiscretize(semi, tspan)

info_callback = InfoCallback(interval=50)

function l_2(system::ParticlePackingSystem, dv_ode, du_ode, v_ode, u_ode, semi, t)
    system.is_boundary && return nothing

    u = TrixiParticles.wrap_u(u_ode, system, semi)
    TrixiParticles.summation_density!(system, semi, u, u_ode, system.density)

    return sqrt(sum((system.density - system.initial_condition.density) .^ 2) /
                nparticles(system))
end

function l_inf(system::ParticlePackingSystem, dv_ode, du_ode, v_ode, u_ode, semi, t)
    system.is_boundary && return nothing

    u = TrixiParticles.wrap_u(u_ode, system, semi)
    TrixiParticles.summation_density!(system, semi, u, u_ode, system.density)

    return maximum(abs.(system.density - system.initial_condition.density))
end

write_error = PostprocessCallback(; interval=10, l2=l_2, linf=l_inf, write_file_interval=10,
                                  output_directory=output_directory)

callbacks = CallbackSet(UpdateCallback(), info_callback, write_error)

sol = solve(ode, RDPK3SpFSAL35();
            abstol=abstol, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
            reltol=reltol, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
            save_everystep=false, maxiters=maxiters_, callback=callbacks)

packed_ic = InitialCondition(sol, packing_system, semi)
packed_ic.density .= packing_system.density
packed_boundary_ic = InitialCondition(sol, boundary_system, semi)

trixi2vtk(packed_boundary_ic, output_directory=output_directory,
          filename="packed_aorta_boundary",
          pressure=boundary_system.signed_distances) # misuse the pressure for `signed_distance_field`

# Extract the other `InitialCondition`s
geometries_shifted = shift_planar_geometry.(geometries[2:end], transition_length)
boundaries_extruded = extrude_geometry.(geometries_shifted, extrusion_length)
for (i, geometry) in enumerate(boundaries_extruded)
    ic = intersect(packed_ic, geometry)

    trixi2vtk(ic, output_directory=output_directory,
              filename="packed_" * geometry_boundaries[i])
end

geometries_blood_transition = extrude_geometry.(geometries[2:end], transition_length,
                                                omit_bottom_face=true)

geometry_blood_domain = union(first(geometries), geometries_blood_transition...)

trixi2vtk(intersect(packed_ic, geometry_blood_domain), output_directory=output_directory,
          filename="packed_aorta")
