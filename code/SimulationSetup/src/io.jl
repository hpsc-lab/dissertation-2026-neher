"""
Directory and path management for SimulationSetup.

Centralizes all file path handling for data, output, and figures.
"""

# Return path to the figures directory (relative to package repository root)
function fig_dir()::String
    return pkgdir(SimulationSetup, "..", "..", "figures")
end

"""
Return path to the data directory.

Supports three fallback options in order:
1. If configured via `set_config!(data_dir=...)`, use that
2. Standard location relative to package: `../data`

Returns the configured or default data directory.
"""
function data_dir()::String
    config = get_config()
    !isnothing(config.data_dir) && return config.data_dir
    return pkgdir(SimulationSetup, "..", "..", "data")
end

"""Return an output directory for results.

If a `set_out_dir()` hook is defined in Main and no `result_variant` is given,
that hook is used. Otherwise a path under the package 'out' folder is returned.
"""
function out_dir(; result_variant=nothing)
    if (isdefined(Main, :set_out_dir) && isnothing(result_variant))
        return Main.set_out_dir()
    end

    run_id_ = isempty(current_run_id()) ? "" : "_$(current_run_id())"

    output_directory = pkgdir(SimulationSetup, "..", "..", "out" * run_id_)

    (isnothing(result_variant)) && return output_directory

    if result_variant == :vulcan || result_variant == :hunter
        output_directory = pkgdir(SimulationSetup, "..", "..", "..",
                                  "out" * "_$result_variant")
        return joinpath(output_directory, "out_v$(current_version())")
    end

    output_directory = pkgdir(SimulationSetup, "..", "..", "out")

    return joinpath(output_directory, "out_$(result_variant)")
end

# Create and return output directory
function ensure_out_dir(; result_variant=nothing)::String
    dir = out_dir(; result_variant=result_variant)
    mkpath(dir)
    return dir
end

# Create and return data directory
function ensure_data_dir()::String
    dir = data_dir()
    mkpath(dir)
    return dir
end

# Create and return figures directory
function ensure_fig_dir()::String
    dir = fig_dir()
    mkpath(dir)
    return dir
end
