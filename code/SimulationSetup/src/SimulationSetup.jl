module SimulationSetup

using Reexport: @reexport
using Suppressor: @suppress

using Dates
using TOML

using LinearAlgebra: normalize, norm, dot, cross, tr, I
using CSV: read
using DataFrames: DataFrames, DataFrame
@reexport using StaticArrays: SVector, SMatrix

using TrixiParticles
using PointNeighbors
@reexport using TrixiParticles: nfaces, TriangleMesh, load_geometry, summary_header,
                                summary_line, summary_footer, @autoinfiltrate,
                                find_too_close_particles, AbstractFluidSystem,
                                nparticles, ndims, wrap_v, wrap_u,
                                current_coordinates, foreach_system, foreach_point_neighbor,
                                hydrodynamic_mass, current_density, viscous_velocity,
                                smoothing_kernel_grad, each_integrated_particle, @threaded,
                                smoothing_kernel, get_neighborhood_search, compact_support,
                                ideal_neighbor_count, default_backend, current_velocity

# Core modules (in dependency order)
include("config.jl")
include("versioning.jl")
include("run_metadata.jl")
include("io.jl")
include("logging.jl")
include("util.jl")

# Preprocessing modules
include("preprocessing/geometric_operations.jl")
include("preprocessing/geometry_io.jl")
include("preprocessing/hemodynamic_parameters.jl")
include("preprocessing/simulation_parameters.jl")

# Postprocessing modules
include("postprocessing/boundary_reporters.jl")
include("postprocessing/wall_shear_stress.jl")
include("postprocessing/turbulence_model.jl")
include("postprocessing/concatenate_results.jl")
include("postprocessing/prepare_plot_data.jl")
include("postprocessing/analyze_slices.jl")

# Export configuration API
export set_config!, get_config, reset_config!

# Export versioning API
export current_version, infer_scenario, current_scenario, current_code_version,
       current_run_id, git_command_output, git_is_dirty, git_branch, git_commit,
       initialize_code_version!

# Export run metadata API
export RunContext, create_run_context, set_run_context!, load_run_context,
       write_run_metadata, timestamp_utc

# Export I/O API
export fig_dir, data_dir, out_dir, ensure_out_dir, ensure_data_dir, ensure_fig_dir

# Export logging API
export print_startup_message, print_restart_message

# Export utilities
export coords_eltype, dict_to_sorted_vector, latest_restart_iter, latest_simulation_iter

# Export configuration components
export ncycles, transition_length
export cm_to_m, cm2_to_m2, m2_to_cm2, ml_to_m3, m3_to_ml, m_to_mm, Pa_to_mmHg, mmHg_to_Pa
export boundary_ids, boundary_names
export subject_parameters, get_parameter_table, extract_subject_parameters

# Export preprocessing
export LumpedParameters, SimulationParameters
export baeumler_flow_ratios, realistic_flow_ratios, unrealistic_flow_ratios,
       adapted_flow_ratios
export shift_planar_geometry, surface_area, project_points_to_plane
export save_stl

# Export postprocessing - boundary reporters
export pp_functions

# Export postprocessing - wall shear stress
export compute_wss_timeseries, wall_shear_stress, SPSTurbulenceModelDalrymple

# Export postprocessing - concatenation
export concatenate_time_series, concatenate_csv

# Export postprocessing - plot data
export prepare_plot_data
export precalculate_pressure

# Export postprocessing - slice analysis
export interpolate_fluid_properties, interpolate_structure_properties, write_slices_to_vtk

end # module SimulationSetup
