# Calculate the inlet velocity using Fourier series approximation,
# Slightly adapted from Zhang et al. (2025), Eq. 21.
# Approximates the continuous cardiac cycle of the waveform for aorta flow.
# V_inlet(t) = a₀ + Σ(n=1 to 8) aₙ cos(nωt) + Σ(n=1 to 8) bₙ sin(nωt)
function velocity_inlet_fourier(t)
    # Fourier coefficients from empirical data (Table III)
    a_0 = 0.308

    # Eq. 21: V_inlet(t) = a₀ + Σ aₙ cos(nωt) + Σ bₙ sin(nωt)
    v_inlet = a_0

    # Cosine and sine terms with coefficients directly embedded
    v_inlet += (-0.1812) * cos(1 * omega * t) + (-0.07725) * sin(1 * omega * t)
    v_inlet += 0.1276 * cos(2 * omega * t) + 0.01466 * sin(2 * omega * t)
    v_inlet += (-0.08981) * cos(3 * omega * t) + 0.04295 * sin(3 * omega * t)
    v_inlet += 0.04347 * cos(4 * omega * t) + (-0.06679) * sin(4 * omega * t)
    v_inlet += (-0.05412) * cos(5 * omega * t) + 0.05679 * sin(5 * omega * t)
    v_inlet += 0.02642 * cos(6 * omega * t) + (-0.01878) * sin(6 * omega * t)
    v_inlet += 0.008946 * cos(7 * omega * t) + 0.01869 * sin(7 * omega * t)
    v_inlet += (-0.009005) * cos(8 * omega * t) + (-0.01888) * sin(8 * omega * t)

    return v_inlet
end

function velocity_inlet_pulsatile_sin(t)
    amplitude = 1
    frequency = 3 / T

    t_periodic = mod(t, T)
    t_periodic > 1 / frequency && return 0

    return amplitude * sin(pi * frequency * t_periodic)^4
end
