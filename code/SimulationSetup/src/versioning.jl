"""
Version and Git integration for SimulationSetup.

Handles version tracking, scenario inference, and Git state detection.
"""

# Return the current package version from config
function current_version()::VersionNumber
    return GLOBAL_CONFIG.version
end

# Return the current scenario from config
function current_scenario()::Symbol
    return GLOBAL_CONFIG.scenario
end

# Return the current code version from config
function current_code_version()::String
    return GLOBAL_CONFIG.code_version
end

# Return the current run ID from config
function current_run_id()::String
    return GLOBAL_CONFIG.run_id
end

"""
Infer scenario from version number.

Version ranges:
- v"1.0.0" to v"1.0.35": :normotensive
- v"1.0.36": :exercise
- v"1.0.37"+: :hypertensive
"""
function infer_scenario(version::VersionNumber)::Symbol
    version <= v"1.0.35" && return :normotensive
    version == v"1.0.36" && return :exercise
    return :hypertensive
end

# Execute a git command in the package repository and return output
function git_command_output(args...)
    repo_root = pkgdir(SimulationSetup, "..")
    cmd = `git -C $repo_root $(args...)`
    output = IOBuffer()
    process = run(pipeline(ignorestatus(cmd), stdout=output, stderr=devnull))
    success(process) || return nothing

    result = strip(String(take!(output)))
    return isempty(result) ? nothing : result
end

# Check if the Git repository has uncommitted changes
function git_is_dirty()::Bool
    repo_root = pkgdir(SimulationSetup, "..")
    cmd = `git -C $repo_root diff --quiet --ignore-submodules --exit-code`
    process = run(ignorestatus(cmd))
    return !success(process)
end

# Get the current Git branch name
function git_branch()::String
    return something(git_command_output("rev-parse", "--abbrev-ref", "HEAD"), "unknown")
end

# Get the current Git commit hash (short form)
function git_commit()::String
    return something(git_command_output("rev-parse", "--short", "HEAD"), "unknown")
end

"""
Initialize code version from multiple sources in order of preference:
1. TRIXI_HEMODYNAMICS_CODE_VERSION environment variable
2. Current Git commit hash
3. Fallback to "unknown"
"""
function initialize_code_version!()
    code_version_env = get(ENV, "TRIXI_HEMODYNAMICS_CODE_VERSION", "")
    if !isempty(code_version_env)
        set_config!(code_version=code_version_env)
        return code_version_env
    end

    commit = git_commit()
    if commit != "unknown"
        set_config!(code_version=commit)
        return commit
    end

    set_config!(code_version="unknown")
    return "unknown"
end
