"""Boundary reporter functions and postprocessing callback assembly.

This file provides a set of reporter functions used by the simulation framework
for measuring and reporting inlet/outlet quantities during solver execution.
These callbacks are assembled via `pp_functions()` and passed into the solver.
"""

inflow_velocity(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function inflow_velocity(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode, semi, t)
    return norm(system.cache.velocity_reference_values[1](0, t))
end

Q_inlet(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function Q_inlet(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode, semi, t)
    return system.cache.boundary_zones_flow_rate[1][]
end

p_outlet_thoracic(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function p_outlet_thoracic(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                           semi, t)
    return system.cache.pressure_reference_values[2].pressure[]
end

Q_outlet_thoracic(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function Q_outlet_thoracic(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                           semi, t)
    return system.cache.boundary_zones_flow_rate[2][]
end

p_outlet_left_common(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function p_outlet_left_common(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                              semi, t)
    return system.cache.pressure_reference_values[3].pressure[]
end

Q_outlet_left_common(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function Q_outlet_left_common(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                              semi, t)
    return system.cache.boundary_zones_flow_rate[3][]
end

p_outlet_left_subclavian(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function p_outlet_left_subclavian(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                                  semi, t)
    return system.cache.pressure_reference_values[4].pressure[]
end

Q_outlet_left_subclavian(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function Q_outlet_left_subclavian(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                                  semi, t)
    return system.cache.boundary_zones_flow_rate[4][]
end

p_outlet_brachiocephalic(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function p_outlet_brachiocephalic(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                                  semi, t)
    return system.cache.pressure_reference_values[5].pressure[]
end

Q_outlet_brachiocephalic(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function Q_outlet_brachiocephalic(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                                  semi, t)
    return system.cache.boundary_zones_flow_rate[5][]
end

p_outlet_right_subclavian(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function p_outlet_right_subclavian(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                                   semi, t)
    return system.cache.pressure_reference_values[5].pressure[]
end

Q_outlet_right_subclavian(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function Q_outlet_right_subclavian(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                                   semi, t)
    return system.cache.boundary_zones_flow_rate[5][]
end

p_outlet_right_common(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function p_outlet_right_common(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                               semi, t)
    return system.cache.pressure_reference_values[6].pressure[]
end

Q_outlet_right_common(system, dv_ode, du_ode, v_ode, u_ode, semi, t) = nothing
function Q_outlet_right_common(system::OpenBoundarySystem, dv_ode, du_ode, v_ode, u_ode,
                               semi, t)
    return system.cache.boundary_zones_flow_rate[6][]
end

function total_volume(system, dv_ode, du_ode, v_ode, u_ode, semi, t)
    v = TrixiParticles.wrap_v(v_ode, system, semi)

    density = view(TrixiParticles.current_density(v, system),
                   TrixiParticles.each_active_particle(system))
    mass = view(system.mass, TrixiParticles.each_active_particle(system))

    return mapreduce(+, density, mass) do rho_i, m_i
        return m_i / rho_i
    end
end

function total_mass_(system, dv_ode, du_ode, v_ode, u_ode, semi, t)
    mass = view(system.mass, TrixiParticles.each_active_particle(system))
    return sum(mass)
end

"""Return a NamedTuple of postprocessing callbacks for the solver.

Depending on the outlet configuration the returned NamedTuple contains functions
for inlet velocity, outlet pressures/flows and global quantities (mass, volume).
These callbacks are expected by the time-stepping/monitoring utilities.
"""
function pp_functions(boundary_dict)
    if haskey(boundary_dict, "brachiocephalic")
        return (v_in=inflow_velocity, Q_inlet=Q_inlet,
                p_outlet_thoracic=p_outlet_thoracic,
                Q_outlet_thoracic=Q_outlet_thoracic,
                p_outlet_left_common=p_outlet_left_common,
                Q_outlet_left_common=Q_outlet_left_common,
                p_outlet_left_subclavian=p_outlet_left_subclavian,
                Q_outlet_left_subclavian=Q_outlet_left_subclavian,
                p_outlet_brachiocephalic=p_outlet_brachiocephalic,
                Q_outlet_brachiocephalic=Q_outlet_brachiocephalic,
                total_mass=total_mass_, total_volume=total_volume)
    elseif haskey(boundary_dict, "right_subclavian") &&
           haskey(boundary_dict, "right_common")
        return (v_in=inflow_velocity, Q_inlet=Q_inlet,
                p_outlet_thoracic=p_outlet_thoracic,
                Q_outlet_thoracic=Q_outlet_thoracic,
                p_outlet_left_common=p_outlet_left_common,
                Q_outlet_left_common=Q_outlet_left_common,
                p_outlet_left_subclavian=p_outlet_left_subclavian,
                Q_outlet_left_subclavian=Q_outlet_left_subclavian,
                p_outlet_right_subclavian=p_outlet_right_subclavian,
                Q_outlet_right_subclavian=Q_outlet_right_subclavian,
                p_outlet_right_common=p_outlet_right_common,
                Q_outlet_right_common=Q_outlet_right_common,
                total_mass=total_mass_, total_volume=total_volume)
    else
        error("unkown configuration")
    end
end
