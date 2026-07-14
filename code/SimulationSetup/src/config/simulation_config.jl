"""
Simulation configuration module.

Organizes and loads all configuration-related modules for cardiovascular
simulation setup, including boundaries, unit conversions, subject data loading,
and simulation constants.
"""

# Core config components (order matters for dependencies)
include("simulation_constants.jl")
include("unit_conversion.jl")
include("boundaries.jl")
include("subject_loader.jl")
