using SimulationSetup
using TrixiParticles

v_x_interpolated(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function v_x_interpolated(system::TrixiParticles.AbstractFluidSystem{3},
                          dv_ode, du_ode, v_ode, u_ode, semi, t)
    start_point = [flow_length / 2, -pipe_radius, 0.0]
    end_point = [flow_length / 2, pipe_radius, 0.0]

    values = interpolate_line(start_point, end_point, 100, semi, system, v_ode, u_ode;
                              cut_off_bnd=true, clip_negative_pressure=false,
                              include_wall_velocity=true)

    return values.velocity[1, :]
end

particle_spacing_factor = 30
output_directory = joinpath(out_dir(), "validation",
                            "open_boundaries", "pulsatile_channel_flow_3d")
filename = "result_vx" * "_dp_$(particle_spacing_factor)"
pp_callback = PostprocessCallback(; dt=0.01, output_directory=output_directory,
                                  v_x=v_x_interpolated, filename=filename,
                                  write_csv=true, write_file_interval=1)

trixi_include(@__MODULE__,
              pkgdir(SimulationSetup, "..", "scripts", "validation", "open_boundaries",
                     "pulsatile_channel_flow_3d.jl"),
              tspan=(0.0, 12.6), extra_callback=pp_callback,
              particle_spacing_factor=particle_spacing_factor)
