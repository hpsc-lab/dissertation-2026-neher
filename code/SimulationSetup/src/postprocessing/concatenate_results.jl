"""Concatenate VTU time-series produced by restart cycles into a single
consistent time series.

This function scans restart_* folders, maps restart-local iterators to a global
iteration grid and writes VTU files with consistent iteration numbers.
"""
function concatenate_time_series(results_dir, semi; tspan=(0.0, 0.75), dt=0.01, FSI=false)
    last_iter_restart = latest_restart_iter(results_dir)

    # No restart files found - nothing to concatenate
    last_iter_restart == 0 && return results_dir

    # Generate time vector for the complete simulation
    times = first(tspan):dt:last(tspan)
    t0 = first(tspan)
    n_times = length(times)
    written_iters = Set{Int}()

    # Map a loaded time to a global iteration index in a tolerance-safe way
    function time_to_iter(t_loaded; atol=1e-6)
        idx = round(Int, (t_loaded - t0) / dt) + 1
        idx = clamp(idx, 1, n_times)
        if isapprox(times[idx], t_loaded; atol=atol)
            return idx - 1
        end
        return nothing
    end

    # Process each restart iteration
    for restart_iteration in 1:last_iter_restart
        @info "Processing restart iteration $restart_iteration/$last_iter_restart"

        # Build file prefix for current restart iteration
        prefix = restart_iteration > 0 ? "restart_$(restart_iteration)_" : ""
        last_iter_simulation = latest_simulation_iter(results_dir, prefix * "fluid_1")

        # Process each simulation iteration within this restart
        for simulation_iteration in 0:last_iter_simulation
            # Build file paths for restart files
            suffix = "_$simulation_iteration.vtu"
            file_fluid = joinpath(results_dir, prefix * "fluid_1" * suffix)
            file_open_boundary = joinpath(results_dir, prefix * "open_boundary_1" * suffix)
            boundary_name = FSI ? "structure_1" : "boundary_1"
            file_boundary = joinpath(results_dir, prefix * boundary_name * suffix)

            # Load restart state
            ode_dummy = semidiscretize(semi, tspan;
                                       restart_with=(file_fluid, file_open_boundary,
                                                     file_boundary))

            # Derive global iteration index from the loaded time
            t_loaded = first(ode_dummy.tspan)
            iter = time_to_iter(t_loaded; atol=1e-6)

            if iter === nothing
                @warn "  ⊳ Skipped file $suffix - time $t_loaded not on expected grid"
                continue
            end

            if iter in written_iters
                @info "  ⊳ Skipped iter $iter (t = $t_loaded s) - already written"
                continue
            end

            @info "  ⊲ Processing iter $iter (t = $t_loaded s)"

            # Verify that loaded time matches expected time
            if isapprox(times[iter + 1], t_loaded; atol=1e-6)
                @info "    ↳ Writing VTU files for iter $iter ($t_loaded sec.)..."
                # Write VTU files with correct global iteration number
                # Use the top-level writer so system quantities (e.g., pressure) are updated.
                trixi2vtk(ode_dummy.u0, semi, t_loaded; output_directory=results_dir,
                          iter=iter)
                push!(written_iters, iter)
                @info "    ✓ Completed iter $iter"
            else
                @warn "  ⊳ Skipped iter $iter - time mismatch (expected $(times[iter + 1]), got $t_loaded)"
            end
        end
    end

    @info "Concatenation complete"
end

"""Concatenate CSV pressure results from restart parts into a single CSV.

The original file is backed up as restart_0_resulting_pressures.csv.
"""
function concatenate_csv(results_dir)
    # Load initial pressure data (restart iteration 0)
    data_sim = [CSV.read(joinpath(results_dir, "resulting_pressures.csv"), DataFrame)]

    # Load all restart pressure data files
    n_iterations = latest_restart_iter(results_dir)
    for iter in 1:n_iterations
        push!(data_sim,
              CSV.read(joinpath(results_dir, "restart_$(iter)_resulting_pressures.csv"),
                       DataFrame))
    end

    # Combine all data and clean up
    combined = vcat(data_sim...)  # Concatenate all DataFrames
    combined.time = round.(combined.time, digits=2)  # Round time to avoid floating point issues
    sort!(combined, :time)  # Sort by time
    unique!(combined, :time; keep=:first)  # Remove duplicate time points, keep first occurrence

    # Backup: Rename original file to restart_0_*
    original_file = joinpath(results_dir, "resulting_pressures.csv")
    renamed_file = joinpath(results_dir, "restart_0_resulting_pressures.csv")
    mv(original_file, renamed_file; force=true)

    # Write combined result with original filename
    CSV.write(joinpath(results_dir, "resulting_pressures.csv"), combined)
end

"""Interactively delete restart files in a results directory.

If `only_vtu=true` only `.vtu` restart files are considered. This function prompts
for confirmation before deleting files.
"""
function delete_restart_files(results_dir; only_vtu=true)
    files = readdir(results_dir)
    restart_files = filter(f -> startswith(f, "restart_"), files)

    if only_vtu
        restart_files = filter(f -> endswith(f, ".vtu"), restart_files)
    end

    if isempty(restart_files)
        @info "no restart files found"
        return
    end

    file_type = only_vtu ? "restart VTU files" : "restart files"
    @warn "Found $(length(restart_files)) $file_type:"
    for file in restart_files
        println("  - $file")
    end

    print("\nDelete all $file_type? (yes/no): ")
    response = readline()

    if lowercase(strip(response)) in ["yes", "y"]
        for file in restart_files
            file_path = joinpath(results_dir, file)
            rm(file_path; force=true)
            @info "deleted restart file: $file"
        end
        @info "finished deleting $file_type"
    else
        @info "deletion cancelled"
    end
end

"""Append a single iteration file to an existing ParaView collection/file.

Used to build PVD collections for visualization time series.
"""
function append_collection(file_vtk, iter)
    file_input = file_vtk * "_iter_$iter.vtu"
    input = vtk2trixi(file_input; element_type=Float64,
                      create_initial_condition=false)
    pvd = TrixiParticles.paraview_collection.(file_vtk; append=true)

    # Create VTK vertex cells for each wall particle
    center_point = [input.coordinates[:, end];;]
    cells = [TrixiParticles.MeshCell(TrixiParticles.VTKCellTypes.VTK_VERTEX, (i,))
             for i in axes(center_point, 2)]
    t = iter * 0.01
    TrixiParticles.vtk_grid(file_input, center_point, cells) do vtk
        # Store basic particle information
        vtk["index"] = TrixiParticles.eachindex(cells)
        vtk["time"] = t
        vtk["ndims"] = 3

        vtk["particle_spacing"] = 0.001

        # Add to ParaView collection
        pvd[t] = vtk
    end

    # Save collection file (only for time series iterations)
    TrixiParticles.vtk_save(pvd)
end
