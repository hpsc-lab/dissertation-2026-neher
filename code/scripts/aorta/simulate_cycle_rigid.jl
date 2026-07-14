using SimulationSetup
using TrixiParticles

# Load variables into scope
trixi_include(@__MODULE__, joinpath(@__DIR__, "setup_rigid.jl");
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

ode_restart.tspan = (zero(tspan[1]), convert(eltype(tspan), T))

output_directory_cycle = joinpath(output_directory, "full_cycle")
saving_cb = SolutionSavingCallback(dt=0.01, overwrite=false,
                                   output_directory=output_directory_cycle)

pp_cb = PostprocessCallback(; dt=0.01, output_directory=output_directory_cycle,
                            filename="resulting_pressures", write_csv=true,
                            write_file_interval=1, pp_functions(boundaries)...)

callbacks = CallbackSet(InfoCallback(interval=50), SortingCallback(; interval=1000),
                        saving_cb, pp_cb, UpdateCallback())

sol_restart = solve(ode_restart, RDPK3SpFSAL35(),
                    abstol=1e-7, # Default abstol is 1e-6 (may need to be tuned to prevent boundary penetration)
                    reltol=1e-4, # Default reltol is 1e-3 (may need to be tuned to prevent boundary penetration)
                    dtmax=1e-2, # Limit stepsize to prevent crashing
                    save_everystep=false, callback=callbacks)
