"""
Logging and message output for simulation startup and restart.

Provides formatted summaries of simulation configuration and state.
"""

"""
Print a summary of simulation startup information.

Generic implementation that works with any parameters object and custom info pairs.
Displays current version, scenario, run ID, and any additional custom information.

# Arguments
- `parameters`: Simulation parameters object to display
- `custom_infos`: Variable number of (key, value) pairs to include in summary
"""
function print_startup_message(parameters; custom_infos...)
    io = stdout
    io_context = IOContext(io,
                           :compact => false,
                           :key_width => 30,
                           :total_width => 100,
                           :indentation_level => 0)

    summary_header(io,
                   "Simulation Startup: " * parameters.subject)
    summary_line(io, "version", current_version())
    summary_line(io, "scenario", current_scenario())
    summary_line(io, "run id", current_run_id())

    # Include any custom info pairs
    for (key, info) in custom_infos
        summary_line(io, string(key), info)
    end

    summary_footer(io)
    println(io, "\n")

    # Display full parameters object
    show(io_context, MIME"text/plain"(), parameters)
    println(io, "\n")
end

"""
Print simulation restart information.

Displays restart status with current configuration and progress information.

# Arguments
- `parameters`: Simulation parameters object
- `latest_iter`: Latest completed iteration number
- `key`: Symbol determining output type (`:finished` or `:notfinished`)
"""
function print_restart_message(parameters, latest_iter; key=:notfinished)
    io = stdout
    io_context = IOContext(io,
                           :compact => false,
                           :key_width => 30,
                           :total_width => 100,
                           :indentation_level => 0)

    if key === :finished
        summary_header(io,
                       "Simulation Completed: " * parameters.subject)
        summary_line(io, "version", current_version())
        summary_line(io, "scenario", current_scenario())
        summary_line(io, "run id", current_run_id())
        summary_line(io, "final iteration", latest_iter)
        summary_footer(io)
    else
        summary_header(io,
                       "Simulation Restarting: " * parameters.subject)
        summary_line(io, "version", current_version())
        summary_line(io, "scenario", current_scenario())
        summary_line(io, "run id", current_run_id())
        summary_line(io, "next iteration", latest_iter + 1)
        summary_footer(io)
    end

    println(io, "\n")
end
