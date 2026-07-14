"""Compute time-resolved wall shear stress (WSS) quantities and write VTK output.

- compute_wss_timeseries: loops over time steps, computes instantaneous WSS and
  accumulates time-averaged and peak-systolic metrics.
"""
function compute_wss_timeseries(semi, turbulence_model_fluid, results_dir, tspan;
                                t_peak=0.4, dt=0.01, FSI=false,
                                output_directory=results_dir)
    # Generate time vector for the complete simulation
    times = tspan[1]:dt:tspan[2]

    @info "Determining sample points on wall surface..."
    wall_particles = surface_particles_from_boundary(semi)
    @info "  ✓ Found $(length(wall_particles)) wall particles"

    # Initialize arrays to accumulate WSS values over time
    WSS_values_vec = fill(zero(SVector{3, eltype(semi.systems[1])}), length(wall_particles))
    WSS_values_mag = zeros(length(wall_particles))

    # Storage for peak systolic wall shear stress (at t = t_peak)
    PSWSS = fill(zero(SVector{3, eltype(semi.systems[1])}), length(wall_particles))

    initial_wall_sampled = InitialCondition[]

    # Process each time step
    for (i, t) in enumerate(times)
        iter = i - 1

        @info "Processing time step $i/$(length(times)) (t = $t s)..."

        # Build file paths for restart files
        file_fluid = joinpath(results_dir, "fluid_1_$iter.vtu")
        file_open_boundary = joinpath(results_dir, "open_boundary_1_$iter.vtu")
        boundary_name = FSI ? "structure_" : "boundary_"
        file_boundary = joinpath(results_dir, boundary_name * "1_$iter.vtu")

        # Load simulation state from restart files
        ode = semidiscretize(semi, (t, tspan[2]);
                             restart_with=(file_fluid, file_open_boundary, file_boundary))

        @info "  ↳ Calculating wall shear stress..."
        wall_sampled, tau_w = wall_shear_stress(ode, turbulence_model_fluid, wall_particles)

        # Accumulate WSS values for time averaging
        WSS_values_vec .+= tau_w
        WSS_values_mag .+= TrixiParticles.norm.(tau_w)

        # Write instantaneous WSS to VTK file
        write_vtk_file(wall_sampled, t; output_directory=output_directory,
                       iter=iter, tau_w=tau_w)

        # Store initial wall configuration for later output
        if iter == 0
            push!(initial_wall_sampled, wall_sampled)
        end

        # Store PSWSS at peak systolic time
        if t == t_peak
            @info "  ★ Peak systolic time reached - storing PSWSS"
            PSWSS .= tau_w
        end

        @info "  ✓ Completed time step $i"
    end

    # Calculate time-averaged WSS quantities
    @info "Computing time-averaged quantities..."
    T = tspan[2] - tspan[1]
    WSS_mean_vec = WSS_values_vec .* dt / T
    WSS_mean_mag = WSS_values_mag .* dt / T

    # Calculate derived hemodynamic indices
    OSI = 0.5 .* (1.0 .- TrixiParticles.norm.(WSS_mean_vec) ./ WSS_mean_mag)  # Oscillatory Shear Index
    ECAP = OSI ./ WSS_mean_mag  # Endothelial Cell Activation Potential

    # Write time-averaged quantities to separate VTK file
    @info "Writing time-averaged results..."
    write_vtk_file(first(initial_wall_sampled), 0; output_directory=output_directory,
                   prefix="time_averaged_", iter=-1, OSI=OSI, ECAP=ECAP, PSWSS=PSWSS,
                   WSS_mean_vec=WSS_mean_vec, WSS_mean_mag=WSS_mean_mag)

    @info "✓ WSS time series computation complete"
    return semi
end

"""Compute wall shear stress at sampling points derived from the wall system.

Returns a NamedTuple with `wall_sampled` (InitialCondition) and `tau_w`
(stress vectors).
"""
function wall_shear_stress(ode, turbulence_model_fluid, wall_particles)
    v_ode, u_ode = ode.u0.x
    semi = ode.p
    system = ode.p.systems[1]

    particle_spacing = first(system.initial_condition.particle_spacing)

    # Update systems and neighborhood search to ensure wall velocity is current
    TrixiParticles.update_systems_and_nhs(v_ode, u_ode, semi, 0.0)

    # Calculate stress tensor for the fluid system
    TrixiParticles.calculate_fluid_stress_tensor!(system, turbulence_model_fluid,
                                                  v_ode, u_ode, semi)

    # Extract wall coordinates and shift them slightly inward for sampling
    wall_system = last(semi.systems)
    u = TrixiParticles.wrap_u(u_ode, wall_system, semi)
    wall_coords = TrixiParticles.current_coordinates(u, wall_system)[:, wall_particles]

    # Extrapolate surface normals from neighboring fluid particles
    wall_normals = extrapolat_normals(wall_coords, turbulence_model_fluid, system,
                                      v_ode, u_ode, semi)

    # Shift sampling points inward along surface normals (half particle spacing)
    for point in axes(wall_coords, 2)
        point_position = TrixiParticles.current_coords(wall_coords, system, point)

        point_position_new = point_position - particle_spacing / 2 * wall_normals[point]
        for dim in 1:ndims(system)
            wall_coords[dim, point] = point_position_new[dim]
        end
    end

    # Create initial condition for wall sampling points
    wall_sampled = InitialCondition(; coordinates=wall_coords,
                                    density=first(system.initial_condition.density),
                                    particle_spacing=particle_spacing)

    # Create turbulence model specifically for wall shear stress calculation
    turbulence_model_wall = SPSTurbulenceModelDalrymple(wall_sampled;
                                                        smallest_length_scale=particle_spacing,
                                                        dynamic_viscosity=turbulence_model_fluid.mu,
                                                        only_wall_shear_stress=true)

    # Calculate wall shear stress at sampling points
    TrixiParticles.calculate_wall_shear_stress!(turbulence_model_wall,
                                                turbulence_model_fluid, system,
                                                v_ode, u_ode, semi)

    return (; wall_sampled, tau_w=turbulence_model_wall.field_variables.stress_vectors)
end

"""Write wall-sampling data and optional custom quantities to VTK files.

Accepts `custom_quantities` as pairs (key=>array) which are added to the VTK.
"""
function write_vtk_file(wall_sampled, t; output_directory="out", iter=-1, prefix="",
                        custom_quantities...)
    mkpath(output_directory)

    # Construct output file paths
    collection_file = joinpath(output_directory, prefix * "wall_shear_stress")
    file = collection_file * (iter >= 0 ? "_$iter" : "")

    # Create or append to ParaView collection file (reset at iter=0)
    pvd = TrixiParticles.paraview_collection(collection_file; append=iter > 0)

    # Create VTK vertex cells for each wall particle
    points = wall_sampled.coordinates
    cells = [TrixiParticles.MeshCell(TrixiParticles.VTKCellTypes.VTK_VERTEX, (i,))
             for i in axes(points, 2)]

    TrixiParticles.vtk_grid(file, points, cells) do vtk
        # Store basic particle information
        vtk["index"] = TrixiParticles.eachparticle(wall_sampled)
        vtk["time"] = t
        vtk["ndims"] = ndims(wall_sampled)

        vtk["particle_spacing"] = fill(first(wall_sampled.particle_spacing),
                                       nparticles(wall_sampled))

        # Add custom quantities (WSS, OSI, ECAP, etc.) if provided
        if !isempty(custom_quantities)
            for (key, quantity) in custom_quantities
                if quantity !== nothing
                    vtk[string(key)] = quantity
                end
            end
        end

        # Add to ParaView collection
        pvd[t] = vtk
    end

    # Save collection file (only for time series iterations)
    iter >= 0 && TrixiParticles.vtk_save(pvd)

    return file
end

"""Select boundary particles that lie close to the fluid to form sampling points.

Returns indices of boundary particles considered part of the surface (within 1.3*h).
"""
function surface_particles_from_boundary(semi)
    system = first(semi.systems)  # Fluid system

    particle_spacing = first(system.initial_condition.particle_spacing)

    fluid_coords = system.initial_condition.coordinates
    boundary_coords = last(semi.systems).initial_condition.coordinates

    # Find boundary particles that are within 1.3*h of fluid particles
    # These are considered "surface" particles for WSS calculation
    return TrixiParticles.find_too_close_particles(boundary_coords, fluid_coords,
                                                   1.3 * particle_spacing)
end

"""Extrapolate outward-pointing normals at wall sampling points using nearby fluid data.

Performs SPH-weighted averaging of cached surface normals from the turbulence model.
"""
function extrapolat_normals(wall_coords, turbulence_model_neighbor, neighbor_system,
                            v_ode, u_ode, semi)
    # Initialize arrays for surface normals and volume weights
    surface_normals = fill(zero(SVector{ndims(neighbor_system), eltype(neighbor_system)}),
                           size(wall_coords, 2))
    volume = zeros(eltype(neighbor_system), size(wall_coords, 2))

    # Extract neighbor system state
    u_neighbor = TrixiParticles.wrap_u(u_ode, neighbor_system, semi)
    v_neighbor = TrixiParticles.wrap_v(v_ode, neighbor_system, semi)
    neighbor_coords = TrixiParticles.current_coordinates(u_neighbor, neighbor_system)

    # Perform neighborhood search and extrapolate normals
    nhs = TrixiParticles.get_neighborhood_search(neighbor_system, semi)
    TrixiParticles.foreach_point_neighbor(wall_coords, neighbor_coords, nhs;
                                          parallelization_backend=semi.parallelization_backend) do point,
                                                                                                   neighbor,
                                                                                                   pos_diff,
                                                                                                   distance
        # Calculate kernel weight for SPH interpolation
        m_b = TrixiParticles.hydrodynamic_mass(neighbor_system, neighbor)
        rho_b = TrixiParticles.current_density(v_neighbor, neighbor_system, neighbor)
        kernel_weight = TrixiParticles.smoothing_kernel(neighbor_system, distance,
                                                        neighbor) * m_b / rho_b

        # Accumulate surface normals from neighboring fluid particles
        surface_normals_neighbor = turbulence_model_neighbor.cache.surface_normals

        surface_normals[point] += surface_normals_neighbor[neighbor]
        volume[point] += kernel_weight
    end

    # Normalize the accumulated surface normals
    for point in axes(wall_coords, 2)
        # Check volume to avoid division by zero
        if volume[point] > eps(eltype(volume))
            # Negate and normalize to get outward-pointing normal
            surface_normals[point] = -normalize(surface_normals[point] / volume[point])
        end
    end

    return surface_normals
end
