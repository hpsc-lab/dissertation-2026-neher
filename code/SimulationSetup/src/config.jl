"""
Global configuration system for SimulationSetup.

Provides the singleton GLOBAL_CONFIG struct and API for managing
version, scenario, run ID, and directory settings across the package.
Organizes and loads all configuration-related submodules.
"""

"""
    SimulationConfig

Global configuration holder for simulation metadata and directory settings.

# Fields
- `version::VersionNumber`: Package/code version
- `scenario::Symbol`: Simulation scenario (:normotensive, :exercise, :hypertensive)
- `run_id::String`: Unique run identifier
- `code_version::String`: Git commit or environment-based version identifier
- `data_dir::Union{String, Nothing}`: Override data directory path (optional)
- `out_dir::Union{String, Nothing}`: Override output directory path (optional)
- `coord_eltype::Type`: Floating-point type for coordinate arrays
"""
mutable struct SimulationConfig
    version::VersionNumber
    scenario::Symbol
    run_id::String
    code_version::String
    data_dir::Union{String, Nothing}
    out_dir::Union{String, Nothing}
    coord_eltype::Type
end

"""Global singleton configuration instance."""
const GLOBAL_CONFIG = SimulationConfig(v"1.0.0",      # version
                                       :normotensive, # scenario
                                       "",            # run_id
                                       "unknown",     # code_version
                                       nothing,       # data_dir
                                       nothing,       # out_dir
                                       Float64)       # coord_eltype

"""
    set_config!(; kwargs...)

Update global configuration with new values.

# Keywords
- `version::VersionNumber`: Update package version
- `scenario::Symbol`: Update scenario identifier
- `run_id::String`: Update run ID
- `code_version::String`: Update code version
- `data_dir::String`: Override data directory
- `out_dir::String`: Override output directory
- `coord_eltype::Type`: Floating-point type for coordinates (Float64 or Float32)

# Example
```julia
set_config!(version=v"1.0.5", scenario=:hypertensive, coord_eltype=Float32)
```
"""
function set_config!(; version=nothing, scenario=nothing, run_id=nothing,
                     code_version=nothing, data_dir=nothing, out_dir=nothing,
                     coord_eltype=nothing)
    !isnothing(version) && (GLOBAL_CONFIG.version = version)
    !isnothing(scenario) && (GLOBAL_CONFIG.scenario = scenario)
    !isnothing(run_id) && (GLOBAL_CONFIG.run_id = run_id)
    !isnothing(code_version) && (GLOBAL_CONFIG.code_version = code_version)
    !isnothing(data_dir) && (GLOBAL_CONFIG.data_dir = data_dir)
    !isnothing(out_dir) && (GLOBAL_CONFIG.out_dir = out_dir)
    !isnothing(coord_eltype) && (GLOBAL_CONFIG.coord_eltype = coord_eltype)
    return GLOBAL_CONFIG
end

"""
    get_config()::SimulationConfig

Retrieve the current global configuration.

Returns:
    SimulationConfig struct with current settings.
"""
function get_config()::SimulationConfig
    return GLOBAL_CONFIG
end

"""
    reset_config!()::SimulationConfig

Reset global configuration to defaults.

Returns the reset GLOBAL_CONFIG.
"""
function reset_config!()::SimulationConfig
    GLOBAL_CONFIG.version = v"1.0.0"
    GLOBAL_CONFIG.scenario = :normotensive
    GLOBAL_CONFIG.run_id = ""
    GLOBAL_CONFIG.code_version = "unknown"
    GLOBAL_CONFIG.data_dir = nothing
    GLOBAL_CONFIG.out_dir = nothing
    GLOBAL_CONFIG.coord_eltype = Float64
    return GLOBAL_CONFIG
end

# Load configuration submodules
include("config/simulation_config.jl")
