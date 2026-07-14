using SimulationSetup
using TrixiParticles
include("velocity_functions.jl")

# ==========================================================================================
# ==== Configuration

# Configure simulation with version and scenario
set_config!(version=v"1.0.35", scenario=:normotensive, coord_eltype=Float64)
initialize_code_version!()

# Scenario-dependent parameters: maps scenario → (T, p_syst, p_diast, stroke_vol, v_peak, v_max)
const SCENARIO_PARAMS = Dict(
    :normotensive => (T=0.75, p_syst=125.0, p_diast=75.0, stroke_volume_factor=1.0,
                      v_peak_factor=1.0, v_max=3.0),
    :exercise => (T=0.4, p_syst=180.0, p_diast=85.0, stroke_volume_factor=1.2,
                  v_peak_factor=2.2, v_max=5.0),
    :hypertensive => (T=0.55, p_syst=200.0, p_diast=120.0, stroke_volume_factor=1.05,
                      v_peak_factor=1.3, v_max=3.0)
)

# Extract parameters for current scenario
scenario = current_scenario()
params = get(SCENARIO_PARAMS, scenario, SCENARIO_PARAMS[:normotensive])

const T = params.T
const omega = 2pi / T

density_blood = 1060.0

# Resolution
particle_spacing = 1e-3

# Patient identifier
subject = "F10"

p_syst = params.p_syst
p_diast = params.p_diast
stroke_volume_factor = params.stroke_volume_factor
v_peak_factor = params.v_peak_factor

flow_rate_correction_factor = current_version() <= v"1.0.21" ? 1.3 : 1.0

param_sim = SimulationParameters(subject; particle_spacing, T, density_blood,
                                 q_prescribed=realistic_flow_ratios,
                                 stroke_volume_factor=stroke_volume_factor,
                                 v_peak_factor=v_peak_factor,
                                 p_syst=p_syst, p_diast=p_diast, L_eff=0.35)
(; boundaries) = param_sim

geometry_names = vcat(["aorta", "aorta_boundary"], collect(keys(boundaries)))
files = joinpath.(data_dir(), "aorta_initial_condition",
                  "v$(current_version().major).$(current_version().minor)",
                  "packed_results_" * subject, "dp_$(particle_spacing)",
                  "packed_" .* geometry_names .* ".vtu")
ic_dict = Dict(zip(geometry_names,
                   getfield.(vtk2trixi.(files; element_type=typeof(particle_spacing),
                                        coordinates_eltype=SimulationSetup.coords_eltype()),
                             :initial_condition)))

# ==========================================================================================
# ==== Experiment Setup
open_boundary_layers = 10

v_max = params.v_max

sound_speed_factor = 10
sound_speed = sound_speed_factor * v_max

for key in keys(ic_dict)
    key == "aorta_boundary" && continue
    ic_dict[key].density .= density_blood
    ic_dict[key].mass .= density_blood * particle_spacing^3
end

factor_buffer_size = 1
n_buffer_particles = factor_buffer_size * nparticles(ic_dict["inflow"])

# ==========================================================================================
# ==== Fluid
smoothing_length = 1.4 * particle_spacing
smoothing_kernel = WendlandC2Kernel{3}()

fluid_density_calculator = ContinuityDensity()

kinematic_viscosity = 0.004 / density_blood

viscosity = ViscosityAdami(nu=kinematic_viscosity)

background_pressure = 10 * sound_speed * density_blood * v_max^2
shifting_technique = TransportVelocityAdami(; background_pressure)

state_equation = StateEquationCole(; sound_speed, reference_density=density_blood,
                                   exponent=1)

density_diffusion = DensityDiffusionMolteniColagrossi(delta=0.1)

fluid_system = WeaklyCompressibleSPHSystem(ic_dict["aorta"];
                                           density_calculator=fluid_density_calculator,
                                           state_equation, smoothing_kernel,
                                           density_diffusion=density_diffusion,
                                           smoothing_length, viscosity=viscosity,
                                           shifting_technique=shifting_technique,
                                           buffer_size=n_buffer_particles)

# ==========================================================================================
# ==== Open boundary
boundary_zones = Dict{String, BoundaryZone}()

# Inflow
face_in, face_normal_in = planar_geometry_to_face(boundaries["inflow"].geometry)
const velocity_direction = copy(-face_normal_in)
# The Fourier-based inlet waveform `velocity_inlet_fourier(t)` reaches a peak
# near 1.0 m/s. To impose the target flow rate, scale the velocity waveform by `v_peak`.
const v_peak = param_sim.subject_parameters.v_peak * flow_rate_correction_factor
reference_velocity = (pos, t) -> (v_peak * velocity_direction * velocity_inlet_fourier(t))
inflow = BoundaryZone(; boundary_face=face_in, face_normal=(-face_normal_in),
                      open_boundary_layers, particle_spacing,
                      sample_points=boundaries["inflow"].sample_points,
                      initial_condition=ic_dict["inflow"], density=density_blood,
                      reference_velocity=reference_velocity,
                      boundary_type=BidirectionalFlow())
boundary_zones["inflow"] = inflow

# Outflows
for key in keys(boundaries)
    key == "inflow" && continue
    face_out, face_normal_out = planar_geometry_to_face(boundaries[key].geometry)
    outflow = BoundaryZone(; boundary_face=face_out, face_normal=(-face_normal_out),
                           open_boundary_layers, particle_spacing,
                           sample_points=boundaries[key].sample_points,
                           reference_pressure=boundaries[key].pressure_model,
                           initial_condition=ic_dict[key], density=density_blood,
                           boundary_type=BidirectionalFlow())
    boundary_zones[key] = outflow
end

# We need to pass a sorted vector of `BoundaryZones`s, since the post process functions
# rely on the index of the boundary zone.
boundary_zones_sorted = SimulationSetup.dict_to_sorted_vector(boundary_zones)

open_boundary = OpenBoundarySystem(boundary_zones_sorted...; fluid_system,
                                   boundary_model=BoundaryModelDynamicalPressureZhang(),
                                   density_diffusion=DensityDiffusionAntuono(delta=0.1),
                                   calculate_flow_rate=true,
                                   buffer_size=n_buffer_particles)

# Check configuration and print pressure model values
SimulationSetup.check_pressure_models!(open_boundary, boundaries)
