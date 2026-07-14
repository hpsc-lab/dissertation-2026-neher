# Model pressure dynamics in healthy, stenosed, and arteriosclerotic vessels under pulsatile flow
# Generates three pressure solutions (sol, sol_stenosis, sol_arteriosclerosis) for comparison
# Results are included and plotted in plots/pressure_model/vessel_model.jl

using OrdinaryDiffEqLowOrderRK

include("ode_functions.jl")

T = 1.0
tspan = (0.0, 20T)
dt = T / 1000

vessel_length = 80.0
vessel_diameter = 3.0
density_blood = 1.06

p0 = 0.0

function pulsatile_velocity_sin(t)
    amplitude = 424.0
    frequency = 2 / T

    t_periodic = mod(t, T)
    t_periodic > 1 / frequency && return 0.0

    return amplitude * sin(pi * frequency * t_periodic)^4
end

params_natural = (R_1=0.05, R_2=0.95, C=1.05,
                  q_func=pulsatile_velocity_sin, dt=dt)

# Increase of the resistance of the non-elastic vessels (after outflow) to model narrowed vessels.
params_stenosis = (R_1=0.05 * 1.75, R_2=0.95, C=1.05,
                   q_func=pulsatile_velocity_sin, dt=dt)

# Decrease compliance to model stiffening of the elastic vessels.
params_arteriosclerosis = (R_1=0.05, R_2=0.95, C=1.05 * 0.5,
                           q_func=pulsatile_velocity_sin, dt=dt)

sol = solve(ODEProblem(pressure_RCR_ode!, [p0], tspan, params_natural), Euler(),
            dt=dt, adaptive=false)

sol_stenosis = solve(ODEProblem(pressure_RCR_ode!, [p0], tspan, params_stenosis),
                     Euler(), dt=dt, adaptive=false)

sol_arteriosclerosis = solve(ODEProblem(pressure_RCR_ode!, [p0], tspan,
                                        params_arteriosclerosis), Euler(), dt=dt,
                             adaptive=false)

q_mean = 0.5 *
         sum(pulsatile_velocity_sin.((19T):dt:(20T)) +
             pulsatile_velocity_sin.((19T - dt):dt:(20T - dt))) * dt / T

p_mean = 0.5 * sum(sol((19T):dt:(20T)) + sol((19T - dt):dt:(20T - dt))) * dt / T
p_mean_stenosis = 0.5 *
                  sum(sol_stenosis((19T):dt:(20T)) + sol_stenosis((19T - dt):dt:(20T - dt))) *
                  dt / T
p_mean_arteriosclerosis = 0.5 *
                          sum(sol_arteriosclerosis((19T):dt:(20T)) +
                              sol_arteriosclerosis((19T - dt):dt:(20T - dt))) * dt / T

R = p_mean / q_mean

A = vessel_diameter^2 * pi / 4
wave_speed = 424.0
k = 0.0075 # (g/(cm⁴·s)) --> (mmHg·s/mL)
Z = wave_speed * density_blood / A * k

t0 = 7.75
t1 = 8.0
p0 = first(sol(t0))
p1 = first(sol(t1))
C = (t1 - t0) / (R * log(p0 / p1))
