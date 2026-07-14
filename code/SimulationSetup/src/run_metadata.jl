"""
Run metadata and context management for simulation tracking.

Provides RunContext struct and metadata persistence to TOML files.
"""

"""
Context information for a simulation run.

# Fields
- `scenario::Symbol`: Simulation scenario (e.g., :normotensive)
- `code_version::String`: Code/commit version identifier
- `run_id::String`: Unique run identifier
- `created_at::String`: ISO 8601 timestamp of run creation
"""
struct RunContext
    scenario::Symbol
    code_version::String
    run_id::String
    created_at::String
end

# Return current UTC timestamp in ISO 8601 format (with Z suffix)
function timestamp_utc()::String
    return Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SSZ")
end

"""
Create a RunContext from current configuration.

# Arguments
- `scenario`: Override scenario (uses config if nothing)
- `run_id`: Override run ID (auto-generated if nothing)
- `code_version`: Override code version (uses config if nothing)

Run ID format: `yyyymmddTHHMMSS_<scenario>` if auto-generated.
"""
function create_run_context(; scenario=nothing, run_id=nothing,
                            code_version=nothing)::RunContext
    scenario_ = isnothing(scenario) ? current_scenario() : scenario
    code_version_ = isnothing(code_version) ? current_version() : code_version

    run_id_ = if isnothing(run_id)
        string(Dates.format(Dates.now(Dates.UTC), dateformat"yyyymmddTHHMMSS"),
               "_", scenario_)
    else
        string(run_id)
    end

    return RunContext(Symbol(scenario_), string(code_version_), run_id_, timestamp_utc())
end

"""Apply RunContext to global configuration."""
function set_run_context!(context::RunContext)::RunContext
    set_config!(scenario=context.scenario,
                run_id=context.run_id,
                code_version=context.code_version)
    return context
end

"""
Load RunContext from metadata.toml file.

Returns `nothing` if file doesn't exist or is malformed.
"""
function load_run_context(metadata_file)::Union{RunContext, Nothing}
    !isfile(metadata_file) && return nothing

    metadata = TOML.parsefile(metadata_file)
    !haskey(metadata, "run") && return nothing

    run = metadata["run"]
    haskey(run, "run_id") || return nothing
    haskey(run, "scenario") || return nothing
    haskey(run, "code_version") || return nothing

    created_at = get(run, "created_at", timestamp_utc())
    return RunContext(Symbol(run["scenario"]), string(run["code_version"]),
                      string(run["run_id"]), string(created_at))
end

# Insert value into dictionary if not nothing
function maybe_insert!(dict::Dict{String, Any}, key::String, value)
    isnothing(value) || (dict[key] = value)
    return dict
end

"""
Write or append simulation metadata to metadata.toml file.

Creates metadata.toml in output_directory with run info, simulation parameters, and git state.
On first creation, stores complete context. On subsequent calls, preserves existing context.

# Arguments
- `output_directory`: Target directory for metadata.toml
- `subject`: Subject identifier (optional)
- `particle_spacing`: Particle spacing parameter (optional)
- `model`: Model name/identifier (optional)
- `result_variant`: Result variant tag (optional)

Returns the RunContext written to metadata.
"""
function write_run_metadata(output_directory; subject=nothing, particle_spacing=nothing,
                            model=nothing, result_variant=nothing)::RunContext
    mkpath(output_directory)
    metadata_file = joinpath(output_directory, "metadata.toml")

    context = load_run_context(metadata_file)
    if isnothing(context)
        context = create_run_context()
        metadata = Dict{String, Any}()

        metadata["run"] = Dict{String, Any}(
            "run_id" => context.run_id,
            "created_at" => context.created_at,
            "scenario" => string(context.scenario),
            "version" => string(current_version()),
            "code_version" => context.code_version
        )

        simulation = Dict{String, Any}()
        maybe_insert!(simulation, "subject", subject)
        maybe_insert!(simulation, "particle_spacing", particle_spacing)
        maybe_insert!(simulation, "model", model)
        maybe_insert!(simulation, "result_variant", result_variant)
        metadata["simulation"] = simulation

        metadata["git"] = Dict{String, Any}(
            "branch" => git_branch(),
            "commit" => context.code_version,
            "dirty" => git_is_dirty()
        )

        open(metadata_file, "w") do io
            TOML.print(io, metadata)
        end
    end

    return set_run_context!(context)
end
