"""
Generic utility functions for SimulationSetup.

Low-level helpers for common operations across the package.
"""

# Element type used for coordinate arrays in this project
function coords_eltype()::Type
    return GLOBAL_CONFIG.coord_eltype
end

# Convert dict to sorted vector by geometry_order (skips missing keys)
function dict_to_sorted_vector(dict; geometry_order=boundary_names())
    sorted_pairs = [(name, dict[name]) for name in geometry_order if haskey(dict, name)]
    return [pair[2] for pair in sorted_pairs]
end

# Find latest restart iterator by scanning restart_*.csv files (returns 0 if none)
function latest_restart_iter(dir)
    files = filter(f -> endswith(f, ".csv"), readdir(dir; join=true))
    iters = map(files) do f
        m = match(r"^restart_(\d+)_", basename(f))
        m === nothing ? 0 : parse(Int, m.captures[1])
    end
    return maximum(iters)
end

# Find latest simulation iteration by scanning VTU files with given prefix
function latest_simulation_iter(dir, fileprefix)
    files = filter(f -> endswith(f, ".vtu"), readdir(dir; join=true))
    iters = map(files) do f
        m = match(Regex("^" * fileprefix * "_(\\d+)"), basename(f))
        m === nothing ? 0 : parse(Int, m.captures[1])
    end
    return maximum(iters)
end
