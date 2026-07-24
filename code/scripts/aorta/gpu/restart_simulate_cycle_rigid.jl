using SimulationSetup
using TrixiParticles
using AMDGPU

parallelization_backend = ROCBackend()

# If `MetalBackend`
# coord_eltype = Float32
coord_eltype = Float64

# Load variables into scope
trixi_include_changeprecision(Float32, @__MODULE__,
                              joinpath(@__DIR__, "..", "setup_rigid.jl"),
                              parallelization_backend=parallelization_backend,
                              coord_eltype=coord_eltype,
                              callbacks=nothing, ode=nothing, sol=nothing)

output_directory_cycle = joinpath(output_directory, "full_cycle")

latest_iter_restart = latest_restart_iter(output_directory_cycle)
prefix = latest_iter_restart > 0 ? "restart_$(latest_iter_restart)_" : ""
latest_iter_simulation = latest_simulation_iter(output_directory_cycle, prefix * "fluid_1")
suffix = "_$(latest_iter_simulation)_current.vtu"

restart_file_fluid = joinpath(output_directory_cycle, prefix * "fluid_1" * suffix)
restart_file_open_boundary = joinpath(output_directory_cycle,
                                      prefix * "open_boundary_1" * suffix)
restart_file_boundary = joinpath(output_directory_cycle, prefix * "boundary_1" * suffix)

tspan_ = (zero(tspan[1]), convert(eltype(tspan), T))
ode_restart = semidiscretize(semi, tspan_;
                             restart_with=(restart_file_fluid,
                                           restart_file_open_boundary,
                                           restart_file_boundary))

is_finished = isapprox(ode_restart.tspan[1], tspan[2])
print_restart_message(param_sim, latest_iter_restart,
                      key=is_finished ? :finished : :not_finished)

restart_prefix = "restart_$(latest_iter_restart +1)"
saving_cb = SolutionSavingCallback(dt=0.01f0, prefix=restart_prefix, overwrite=false,
                                   output_directory=output_directory_cycle)

pp_cb = PostprocessCallback(; dt=0.01f0, output_directory=output_directory_cycle,
                            filename=restart_prefix * "_resulting_pressures",
                            write_csv=true, write_file_interval=1,
                            pp_functions(boundaries)...)

callbacks = CallbackSet(InfoCallback(interval=20), SortingCallback(; interval=1000),
                        saving_cb, pp_cb, UpdateCallback())

sol_restart = solve(ode_restart, RDPK3SpFSAL35(),
                    abstol=1.0f-7, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
                    reltol=1.0f-4, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
                    dtmax=1.0f-2, # Limit stepsize to prevent crashing
                    save_everystep=false, callback=callbacks)
