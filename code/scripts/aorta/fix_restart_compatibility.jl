# This script reads VTU files from simulation output and writes missing fields back into them.
# This is essential for the restart method to work properly, as it ensures all required
# data fields are present in the restart files.

using SimulationSetup
using TrixiParticles
import SimulationSetup: CSV, DataFrame

# Configure version and initialize code version tracking
set_config!(version=v"1.0.13")
initialize_code_version!()

iter = 11

# Patient identifier
subject = "F10"
# Resolution
particle_spacing = 0.0005

trixi_include(@__MODULE__, joinpath(@__DIR__, "setup_fluid.jl"), subject=subject,
              flow_rate_correction_factor=1.3, factor_buffer_size=4,
              particle_spacing=particle_spacing)

output_directory = joinpath(out_dir(), "aorta", "$subject", "elastic",
                            "dp_$(particle_spacing)_t_0.002")

keys_Q = [(Symbol(:Q, i), "Q_$(i)") for i in 1:(param_sim.n_outlets + 1)]
data = vtk2trixi(joinpath(output_directory, "restart_$(iter)_open_boundary_1_0.vtu");
                 keys_Q...)
data_sim = CSV.read(joinpath(output_directory, "restart_$(iter)_resulting_pressures.csv"),
                    DataFrame)
p_outlets = Tuple[]
for key in keys(param_sim.boundaries)
    key == "inflow" && continue
    p = first(data_sim[!, "p_outlet_$(key)_open_boundary_1"])
    push!(p_outlets,
          ("boundary_zone_pressure_$(param_sim.boundaries[key].id+1)", p))
end

keys_p_write = [(Symbol(key_), val) for (key_, val) in p_outlets]

key_to_write = Tuple[]
for key in keys(data)
    if key == :coordinates || key == :initial_condition
        continue
    end

    if key == :particle_spacing
        particle_spacing_ = data.particle_spacing * ones(nparticles(data.initial_condition))
        push!(key_to_write, (key, particle_spacing_))
    else
        push!(key_to_write, (key, data[key]))
    end
end

zone_ids = zeros(Int, nparticles(data.initial_condition))
for particle in 1:nparticles(data.initial_condition)
    particle_coords = TrixiParticles.current_coords(data.coordinates, open_boundary,
                                                    particle)

    for (zone_id, boundary_zone) in enumerate(open_boundary.boundary_zones)
        # Check if boundary particle is in the boundary zone
        if TrixiParticles.is_in_boundary_zone(boundary_zone, particle_coords)
            zone_ids[particle] = zone_id
        end
    end
    @assert zone_ids[particle]!=0 "No boundary zone found for particle"
end

trixi2vtk(data.coordinates;
          filename=joinpath(output_directory,
                            "restart_$(iter)_open_boundary_1.vtu"),
          zone_id=zone_ids, vcat(key_to_write, keys_p_write)...)
