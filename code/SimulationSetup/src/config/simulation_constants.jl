"""
Simulation time and configuration constants.

Provides fundamental simulation parameters that are scenario-independent.
"""

"""Return the number of cardiac cycles to simulate."""
 ncycles() = 7

"""Return the boundary transition length based on code version.

For versions > 0.2, uses 5 milliseconds; earlier versions use 0.
"""
function transition_length()
    return current_version() > v"0.2" ? 5e-3 : 0
end
