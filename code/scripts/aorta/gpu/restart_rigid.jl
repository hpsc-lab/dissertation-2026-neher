using SimulationSetup
using TrixiParticles
using Metal

parallelization_backend = MetalBackend()
coord_eltype = parallelization_backend isa MetalBackend ? Float32 : Float64

# Load variables into scope
trixi_include_changeprecision(Float32, @__MODULE__,
                              joinpath(@__DIR__, "..", "setup_rigid.jl");
                              parallelization_backend=parallelization_backend,
                              coord_eltype=coord_eltype,
                              callbacks=nothing, ode=nothing, sol=nothing)

latest_iter = latest_restart_iter(output_directory)
prefix = latest_iter > 0 ? "restart_$(latest_iter)_" : ""
restart_file_fluid = joinpath(output_directory, prefix * "fluid_1_current.vtu")
restart_file_open_boundary = joinpath(output_directory,
                                      prefix * "open_boundary_1_current.vtu")
restart_file_boundary = joinpath(output_directory, prefix * "boundary_1_current.vtu")

ode_restart = semidiscretize(semi, tspan;
                             restart_with=(restart_file_fluid,
                                           restart_file_open_boundary,
                                           restart_file_boundary))

is_finished = isapprox(ode_restart.tspan[1], tspan[2])
print_restart_message(param_sim, latest_iter,
                      key=is_finished ? :finished : :not_finished)

restart_prefix = "restart_$(latest_iter +1)"
saving_cb_restart = SolutionSavingCallback(dt=0.01f0, prefix=restart_prefix, overwrite=true,
                                           output_directory=output_directory)

pp_cb_restart = PostprocessCallback(; dt=0.01f0, output_directory=output_directory,
                                    filename=restart_prefix * "_resulting_pressures",
                                    write_csv=true, write_file_interval=1,
                                    pp_functions(boundaries)...)

callbacks = CallbackSet(InfoCallback(interval=20), SortingCallback(; interval=1000),
                        saving_cb_restart, pp_cb_restart, UpdateCallback())

sol_restart = solve(ode_restart, RDPK3SpFSAL35(),
                    abstol=1.0f-7, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
                    reltol=1.0f-4, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
                    dtmax=1.0f-2, # Limit stepsize to prevent crashing
                    save_everystep=false, callback=callbacks)
