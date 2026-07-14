# Validate Windkessel pressure models (RC and RCR) against experimental data
# Solves pressure ODEs for 2-element (RC) and 3-element (RCR) models
# Results (sol_2_element, sol_3_element) are included and plotted in plots/pressure_model/validation.jl (Fig. 2.6)

using OrdinaryDiffEqLowOrderRK

include("table_data.jl")
include("ode_functions.jl")

T = 0.375
dt = T / 1000
tspan = (0.0, 10T)
p0 = 78.0

params_RC = (R=108.42, C=1.1808e-2, q_func=pulsatile_flow, dt=dt)
params_RCR = (R_2=106.66, R_1=1.7714, C=1.1808e-2, q_func=pulsatile_flow, dt=dt)

sol_3_element = solve(ODEProblem(pressure_RCR_ode!, [p0], tspan, params_RCR), Euler(),
                      dt=dt,
                      adaptive=false)
sol_2_element = solve(ODEProblem(pressure_RC_ode!, [p0], tspan, params_RC), Euler(), dt=dt,
                      adaptive=false)
