"""Assemble plotting-ready time series from simulation outputs.

Extracts inlet/outlet flows and pressures, computes cycle-averaged values and
converts units to ml/min and mmHg for plotting. Returns a NamedTuple with
(dt, times, p, Q, Q_total, p_mean, Q_mean, p_mean_cycle, Q_mean_cycle).
"""
function prepare_plot_data(data_sim, param_sim, times; T=0.75, t_final=last(times))
    dt = times[2] - times[1]
    n_times = length(times)
    boundaries = param_sim.boundaries

    # Inlet flow: stored as Q_in or derived from inlet velocity
    if hasproperty(data_sim, :Q_inlet_open_boundary_1)
        Q_in = -data_sim[!, "Q_inlet_open_boundary_1"]
    else
        Q_in = data_sim[!, "v_in_open_boundary_1"] .* param_sim.subject_parameters.A_in
    end

    Q_outlets = zeros(param_sim.n_outlets, n_times)
    p_outlets = zeros(param_sim.n_outlets, n_times)

    for key in keys(boundaries)
        key == "inflow" && continue
        outlet_id = boundaries[key].id
        Q_outlets[outlet_id, :] .= data_sim[!, "Q_outlet_$(key)_open_boundary_1"][1:n_times]
        p_outlets[outlet_id, :] .= data_sim[!, "p_outlet_$(key)_open_boundary_1"][1:n_times]
    end

    # Mean values for each cycle
    t0 = times[1]
    cycle_edges = collect(t0:T:t_final)
    if length(cycle_edges) < 2
        cycle_edges = [t_final - T, t_final]
    end

    n_cycles = length(cycle_edges) - 1
    Q_mean_outlets = zeros(param_sim.n_outlets, n_cycles)
    p_mean_outlets = zeros(param_sim.n_outlets, n_cycles)
    Q_mean_t = zeros(param_sim.n_outlets, n_times)
    p_mean_t = zeros(param_sim.n_outlets, n_times)
    for cycle_idx in 1:n_cycles
        t_start = cycle_edges[cycle_idx]
        t_end = cycle_edges[cycle_idx + 1]
        t_mask = t_start .<= times .<= t_end
        if count(t_mask) < 2
            continue
        end
        for key in keys(boundaries)
            key == "inflow" && continue
            outlet_id = boundaries[key].id
            Q_ = Q_outlets[outlet_id, t_mask]
            p_ = p_outlets[outlet_id, t_mask]
            Q_cycle = sum(Q_[2:end] + Q_[1:(end - 1)]) * dt / 2T
            p_cycle = sum(p_[2:end] + p_[1:(end - 1)]) * dt / 2T
            Q_mean_outlets[outlet_id, cycle_idx] = Q_cycle
            p_mean_outlets[outlet_id, cycle_idx] = p_cycle
            Q_mean_t[outlet_id, t_mask] .= Q_cycle
            p_mean_t[outlet_id, t_mask] .= p_cycle
        end
    end

    # Convert units: flow from m³/s to ml/min, pressure from Pa to mmHg
    Q_outlets .*= m3_to_ml()
    Q_in .*= m3_to_ml()
    Q_mean_outlets .*= m3_to_ml()
    Q_mean_t .*= m3_to_ml()
    p_outlets .*= Pa_to_mmHg()
    p_mean_outlets .*= Pa_to_mmHg()
    p_mean_t .*= Pa_to_mmHg()

    return (dt=dt, times=times, p=p_outlets, Q=Q_outlets, Q_total=Q_in,
            p_mean=p_mean_t, Q_mean=Q_mean_t,
            p_mean_cycle=p_mean_outlets, Q_mean_cycle=Q_mean_outlets)
end

"""Precompute boundary pressures and flows for Windkessel models.

Solves RCR-type ODEs for each outlet using the provided solver and returns a
collection of time series (p, Q) converted into plotting units.
"""
function precalculate_pressure(param_sim, times; T=0.75, t_final=ncycles() * T, solve_func,
                               ode_problem, time_integrator, velocity_function)
    dt_sim = times[2] - times[1]
    n_times = length(times)
    results = (dt=dt_sim, times=times, p=zeros(param_sim.n_outlets, n_times),
               Q=zeros(param_sim.n_outlets, n_times), Q_total=zeros(n_times))

    function pressure_RCR_ode!(dp, p, params, t)
        (; R_1, R_2, C, q_func, dt) = params
        # RCR outlet model with numerical inflow derivative
        dq_dt = (q_func(t) - q_func(t - dt)) / dt
        dp[1] = -p[1] / (R_2 * C) + R_1 * dq_dt + q_func(t) * (R_2 + R_1) / (R_2 * C)
    end

    tspan_sim = (0, t_final)
    boundaries = param_sim.boundaries
    (; v_peak, A_in) = param_sim.subject_parameters
    for key in keys(boundaries)
        key == "inflow" && continue

        boundary = boundaries[key]
        v_peak_i = (v_peak * A_in * param_sim.lumped_parametes.q[key]) /
                   boundary.cross_sectional_area

        R1 = boundary.pressure_model.characteristic_resistance
        R2 = boundary.pressure_model.peripheral_resistance
        C = boundary.pressure_model.compliance

        params_RCR = (R_2=R2, R_1=R1, C=C,
                      q_func=t -> v_peak_i * velocity_function(t) *
                                  boundary.cross_sectional_area, dt=1e-3)
        sol_RCR = solve_func(ode_problem(pressure_RCR_ode!, [0.0], tspan_sim, params_RCR),
                             time_integrator, dt=1e-3, adaptive=false)

        results.p[boundary.id, :] .= vec(stack(sol_RCR(times)))
        results.Q[boundary.id, :] .= params_RCR.q_func.(times)
    end

    results.Q_total .= v_peak * A_in * velocity_function.(times)

    # Mean values for each cycle
    t0 = times[1]
    cycle_edges = collect(t0:T:t_final)
    if length(cycle_edges) < 2
        cycle_edges = [t_final - T, t_final]
    end

    n_cycles = length(cycle_edges) - 1
    p_mean_cycle = zeros(param_sim.n_outlets, n_cycles)
    Q_mean_cycle = zeros(param_sim.n_outlets, n_cycles)
    p_mean_t = zeros(param_sim.n_outlets, n_times)
    Q_mean_t = zeros(param_sim.n_outlets, n_times)
    for cycle_idx in 1:n_cycles
        t_start = cycle_edges[cycle_idx]
        t_end = cycle_edges[cycle_idx + 1]
        t_mask = t_start .<= times .<= t_end
        if count(t_mask) < 2
            continue
        end
        for key in keys(boundaries)
            key == "inflow" && continue
            outlet_id = boundaries[key].id
            p_ = results.p[outlet_id, t_mask]
            Q_ = results.Q[outlet_id, t_mask]
            p_cycle = sum(p_[2:end] + p_[1:(end - 1)]) * dt_sim / 2T
            Q_cycle = sum(Q_[2:end] + Q_[1:(end - 1)]) * dt_sim / 2T
            p_mean_cycle[outlet_id, cycle_idx] = p_cycle
            Q_mean_cycle[outlet_id, cycle_idx] = Q_cycle
            p_mean_t[outlet_id, t_mask] .= p_cycle
            Q_mean_t[outlet_id, t_mask] .= Q_cycle
        end
    end

    results = merge(results,
                    (p_mean=p_mean_t, Q_mean=Q_mean_t,
                     p_mean_cycle=p_mean_cycle, Q_mean_cycle=Q_mean_cycle))

    # Convert units: flow from m³/s to ml/min, pressure from Pa to mmHg
    results.Q .*= m3_to_ml()
    results.Q_total .*= m3_to_ml()
    results.Q_mean .*= m3_to_ml()
    results.Q_mean_cycle .*= m3_to_ml()
    results.p .*= Pa_to_mmHg()
    results.p_mean .*= Pa_to_mmHg()
    results.p_mean_cycle .*= Pa_to_mmHg()

    return results
end
