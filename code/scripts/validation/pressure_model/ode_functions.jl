# The reference pressure values are computed using an ODE that describes
# the behavior of the RCR Windkessel model. The governing equation is:
#
#   dp/dt + p / (R_1 * C) = R_2 * dq/dt (R_1 + R_2) / (R_1 * C) * q
#
# where
#   - p: pressure
#   - q: time-dependent flow, provided by `pulsatile_flow`
#   - R_1: characteristic resistance
#   - R_2: peripheral resistance
#   - C: compliance
#
# The function `pressure_RCR_ode!` implements this equation for numerical solution:
function pressure_RCR_ode!(dp, p, params, t)
    (; R_1, R_2, C, q_func, dt) = params
    dq_dt = (q_func(t) - q_func(t - dt)) / dt  # numerical derivative of inflow
    dp[1] = -p[1] / (R_2 * C) + R_1 * dq_dt + q_func(t) * (R_2 + R_1) / (R_2 * C)
end

function pressure_RC_ode!(dp, p, params, t)
    (; R, C, q_func) = params
    dp[1] = (1 / C) * q_func(t) - (p[1] / (R * C))
end
