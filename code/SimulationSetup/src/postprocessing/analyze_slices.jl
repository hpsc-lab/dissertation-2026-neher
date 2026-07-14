"""High-level wrapper: interpolate fluid properties for multiple slices.

Accepts an array of slice geometries and returns interpolated fluid (and
optionally structure) results for each slice.
"""
function interpolate_fluid_properties(slices, semi_0; fsi=false, iter, sampling_resolution,
                                      padding=0.001, results_directory, wall_thickness,
                                      interpolate_structure=false)
    @info "Starting fluid slice interpolation" iter fsi sampling_resolution

    # Load properties into scope
    file_fluid = joinpath(results_directory, "fluid_1_$iter.vtu")
    file_open_boundary = joinpath(results_directory, "open_boundary_1_$iter.vtu")
    file_wall = joinpath(results_directory,
                         fsi ? "structure_1_$iter.vtu" : "boundary_1_$iter.vtu")

    ode = semidiscretize(semi_0, (0.0, 100.0);
                         restart_with=(file_fluid, file_open_boundary, file_wall))

    semi = ode.p
    v_ode, u_ode = ode.u0.x

    @info "Updating systems and neighborhood searches"
    TrixiParticles.update_systems_and_nhs(v_ode, u_ode, semi, first(ode.tspan))

    results_vector = NamedTuple[]
    results_vector_structure = NamedTuple[]

    for (i, slice) in enumerate(slices)
        @info "procesing slice $i out of $(length(slices))"
        results = interpolate_fluid_properties(slice, semi_0, v_ode, u_ode, fsi,
                                               sampling_resolution, padding)
        push!(results_vector, results)
    end

    if interpolate_structure
        for (i, slice) in enumerate(slices)
            @info "procesing slice $i out of $(length(slices)) (structure)"
            results = interpolate_structure_properties(slice, semi, v_ode, u_ode,
                                                       wall_thickness, sampling_resolution,
                                                       padding)
            push!(results_vector_structure, results)
        end
    end

    return results_vector, results_vector_structure
end

function interpolate_fluid_properties(slice::TriangleMesh, semi, v_ode, u_ode, fsi,
                                      sampling_resolution, padding)
    @info "Creating sampling points"
    if fsi
        # Create planar face from slice.
        rectangular_face, face_normal = planar_geometry_to_face(slice)

        particle_spacing = first(first(semi.systems).initial_condition.particle_spacing)

        u_fluid = TrixiParticles.wrap_u(u_ode, first(semi.systems), semi)
        coords_fluid = TrixiParticles.current_coordinates(u_fluid, first(semi.systems))

        u_wall = TrixiParticles.wrap_u(u_ode, last(semi.systems), semi)
        coords_wall = TrixiParticles.current_coordinates(u_wall, last(semi.systems))

        # Find particles close to the slice boundary and project them onto the slice plane.
        ids = find_too_close_particles(TrixiParticles.initial_coordinates(last(semi.systems)),
                                       stack(slice.vertices), particle_spacing)

        points = project_points_to_plane(coords_wall[:, ids], face_normal, rectangular_face)

        sample_points = create_sample_points(points, sampling_resolution;
                                             intersection_mask=coords_fluid,
                                             particle_spacing=particle_spacing)

        sampled_plane = create_sample_points(points, sampling_resolution;
                                             padding=padding, intersection_mask=nothing,
                                             particle_spacing=particle_spacing)

    else
        sample_points = create_sample_points(slice, sampling_resolution)
        sampled_plane = create_sample_points(slice, sampling_resolution; clip=false,
                                             padding=padding)
    end

    @info "Running fluid interpolation" n_sample_points=size(sample_points, 2)
    results = interpolate_points(sample_points, semi, first(semi.systems), v_ode, u_ode;
                                 include_wall_velocity=true, cut_off_bnd=true,
                                 clip_negative_pressure=false)

    results_plane = set_plane_values(sampled_plane, sample_points, results,
                                     sampling_resolution, slice)

    circle_coords = sample_points
    center = SVector{3, eltype(circle_coords)}(Main.mean(circle_coords, dims=2))
    radii = zeros(size(circle_coords, 2))

    for point in axes(circle_coords, 2)
        dist = norm(SVector{3, eltype(circle_coords)}(circle_coords[:, point]) -
                    SVector(center))
        radii[point] = dist
    end

    mean_radius = Main.mean(radii)
    max_radius = maximum(radii)
    min_radius = minimum(radii)

    center_coords = hcat(slice.min_corner, slice.max_corner, center)

    return (; results..., sampling_resolution, results_plane, mean_radius, max_radius,
            min_radius, center_coords=center_coords)
end

"""Interpolate structural (wall) properties for a single slice.

Projects wall particles onto the slice plane, creates sampling points within the
slice thickness, and computes per-sample displacements, velocities and stresses.
Returns a NamedTuple with interpolation results and plane values.
"""
function interpolate_structure_properties(slice::TriangleMesh, semi, v_ode, u_ode,
                                          wall_thickness, sampling_resolution, padding)
    particle_spacing = first(first(semi.systems).initial_condition.particle_spacing)

    u_wall = TrixiParticles.wrap_u(u_ode, last(semi.systems), semi)
    coords_wall = TrixiParticles.current_coordinates(u_wall, last(semi.systems))

    # Create planar face from slice.
    rectangular_face, face_normal = planar_geometry_to_face(slice)

    @info "Creating sampling points"
    # Find particles close to the slice boundary and project them onto the slice plane.
    ids = find_too_close_particles(TrixiParticles.initial_coordinates(last(semi.systems)),
                                   stack(slice.vertices), particle_spacing)

    # Expand to particles close to the full slice thickness.
    ids_slice = find_too_close_particles(coords_wall, coords_wall[:, ids], wall_thickness)
    points = project_points_to_plane(coords_wall[:, ids_slice], face_normal,
                                     rectangular_face)

    sample_points = create_sample_points(points, sampling_resolution;
                                         intersection_mask=coords_wall,
                                         particle_spacing=particle_spacing)

    sampled_plane = create_sample_points(points, sampling_resolution;
                                         padding=padding, intersection_mask=nothing,
                                         particle_spacing=particle_spacing)

    @info "Running structure interpolation" n_sample_points=size(sample_points, 2)
    results = interpolate_structure_properties(sample_points, semi, last(semi.systems),
                                               v_ode, u_ode)

    results_plane = set_plane_values(sampled_plane, sample_points, results,
                                     sampling_resolution, slice; fluid=false)

    circle_coords = coords_wall[:, ids]
    center = SVector{3, eltype(circle_coords)}(Main.mean(circle_coords, dims=2))
    radii = zeros(size(circle_coords, 2))

    for point in axes(circle_coords, 2)
        dist = norm(SVector{3, eltype(circle_coords)}(circle_coords[:, point]) -
                    SVector(center))
        radii[point] = dist
    end

    mean_radius = Main.mean(radii)
    max_radius = maximum(radii)
    min_radius = minimum(radii)

    return (; results..., sampling_resolution, mean_radius, max_radius, min_radius, center,
            results_plane)
end

"""Interpolate structural properties at given sample points.

Given sample coordinates and a structural system, returns per-point
neighbor counts, displacements, velocities and von Mises stress values.
"""
function interpolate_structure_properties(sample_points::AbstractArray, semi, system,
                                          v_ode, u_ode)
    (; parallelization_backend) = semi

    n_sample_points = size(sample_points, 2)

    u = TrixiParticles.wrap_u(u_ode, system, semi)
    v = TrixiParticles.wrap_v(v_ode, system, semi)
    system_coords = TrixiParticles.current_coordinates(u, system)

    search_radius = TrixiParticles.compact_support(system, system)
    nhs = GridNeighborhoodSearch{ndims(system)}(; search_radius,
                                                n_points=size(system_coords, 2))
    TrixiParticles.PointNeighbors.initialize!(nhs, sample_points, system_coords)

    shepard_coefficient = zeros(n_sample_points)
    neighbor_count = zeros(Int, n_sample_points)

    velocities = zeros(ndims(system), n_sample_points)
    displacements = copy(velocities)
    von_mises = zeros(n_sample_points)

    TrixiParticles.foreach_point_neighbor(sample_points, system_coords, nhs;
                                          parallelization_backend) do point, neighbor,
                                                                      pos_diff, distance
        neighbor_count[point] += 1
        m_b = TrixiParticles.hydrodynamic_mass(system, neighbor)
        volume_b = m_b / TrixiParticles.current_density(v, system, neighbor)
        W_ab = TrixiParticles.smoothing_kernel(system, distance, neighbor)

        velocity = TrixiParticles.current_velocity(v, system, neighbor)
        displ = TrixiParticles.current_coords(system, neighbor) -
                TrixiParticles.initial_coords(system, neighbor)

        shepard_coefficient[point] += volume_b * W_ab

        for i in axes(velocities, 1)
            velocities[i, point] += velocity[i] * volume_b * W_ab
            displacements[i, point] += displ[i] * volume_b * W_ab
        end

        von_mises[point] += TrixiParticles.von_mises_stress(system, neighbor) * volume_b *
                            W_ab
    end

    for point in axes(sample_points, 2)
        if neighbor_count[point] > 0
            for i in axes(velocities, 1)
                velocities[i, point] /= shepard_coefficient[point]
                displacements[i, point] /= shepard_coefficient[point]
            end

            von_mises[point] /= shepard_coefficient[point]
        else
            for i in axes(velocities, 1)
                velocities[i, point] = NaN
                displacements[i, point] = NaN
            end

            von_mises[point] = NaN
        end
    end

    return (; point_coords=sample_points, neighbor_count=neighbor_count,
            displacement=displacements, velocity=velocities, von_mises_stress=von_mises)
end

function set_plane_values(sampled_plane, sample_points, results, sampling_resolution,
                          slice; fluid=true)
    rectangular_face, face_normal = planar_geometry_to_face(slice)
    n_plane_points = size(sampled_plane, 2)
    ids_valued = find_too_close_particles(sampled_plane, sample_points,
                                          0.75sampling_resolution)

    @assert length(ids_valued)==size(sample_points, 2) "found $(length(ids_valued)) instead of $(size(sample_points, 2))"

    if fluid
        density = NaN .* ones(n_plane_points)
        pressure = NaN .* ones(n_plane_points)
        velocity = NaN .* ones(3, n_plane_points)

        density[ids_valued] = results.density
        pressure[ids_valued] = results.pressure
        velocity[:, ids_valued] = results.velocity

        return (; density, pressure, velocity,
                point_coords=map_points_to_xy_plane(sampled_plane, rectangular_face,
                                                    face_normal))
    else
        von_mises_stress = NaN .* ones(n_plane_points)
        velocity = NaN .* ones(3, n_plane_points)
        displacement = NaN .* ones(3, n_plane_points)

        von_mises_stress[ids_valued] = results.von_mises_stress
        velocity[:, ids_valued] = results.velocity
        displacement[:, ids_valued] = results.displacement

        return (; von_mises_stress, displacement, velocity,
                point_coords=map_points_to_xy_plane(sampled_plane, rectangular_face,
                                                    face_normal))
    end
end

function map_points_to_xy_plane(sampled_plane, rectangular_face, normal)
    n_plane_points = size(sampled_plane, 2)

    # Build a strict orthonormal basis from the rectangular face.
    T = eltype(sampled_plane)
    origin = SVector{3, T}(rectangular_face[1])
    edge1 = SVector{3, T}(rectangular_face[2]) - origin
    edge2 = SVector{3, T}(rectangular_face[3]) - origin

    # normal = normalize(cross(edge1, edge2))
    e1 = normalize(edge1)
    e2 = normalize(cross(normal, e1))

    # Project each point onto the plane and map it into the 2D basis.
    coords_2d = zeros(2, n_plane_points)
    for i in 1:n_plane_points
        point = SVector{3, T}(sampled_plane[:, i])
        relative = point - origin
        relative_on_plane = relative - dot(relative, normal) * normal
        coords_2d[1, i] = dot(relative_on_plane, e1)
        coords_2d[2, i] = dot(relative_on_plane, e2)
    end

    # Center the 2D coordinates.
    mean_x = Main.mean(coords_2d[1, :])
    mean_y = Main.mean(coords_2d[2, :])
    centered = coords_2d .- [mean_x; mean_y]

    # Align to principal direction for a consistent in-plane orientation.
    rotated_coords = centered
    if n_plane_points > 1
        cov_xx = sum(centered[1, :] .^ 2) / (n_plane_points - 1)
        cov_yy = sum(centered[2, :] .^ 2) / (n_plane_points - 1)
        cov_xy = sum(centered[1, :] .* centered[2, :]) / (n_plane_points - 1)

        angle = 0.5 * atan(2 * cov_xy, cov_xx - cov_yy)
        cos_a = cos(-angle)
        sin_a = sin(-angle)
        rot_matrix = TrixiParticles.SMatrix{2, 2}(cos_a, sin_a, -sin_a, cos_a)
        rotated_coords = rot_matrix * centered
    end

    # Shift to start from origin.
    min_x = minimum(rotated_coords[1, :])
    min_y = minimum(rotated_coords[2, :])

    # Return as 3D coordinates with z = 0.
    mapped_coords = zeros(3, n_plane_points)
    mapped_coords[1, :] .= rotated_coords[1, :] .- min_x
    mapped_coords[2, :] .= rotated_coords[2, :] .- min_y

    return mapped_coords
end

"""Write VTK files for a set of anatomical slices across simulation iterations.

Iterates over configured slice geometries, runs interpolation (fluid and optional
structure) and writes VTK/InitialCondition files for visualization.
"""
function write_slices_to_vtk(subject, semi, results_directory, sampling_resolution;
                             fsi=false, interpolate_structure=false, wall_thickness,
                             start_iter=0)
    @info "Starting VTK slice writing" fsi sampling_resolution

    n_sim_iters = if current_version() == v"1.0.35"
        75
    elseif current_version() == v"1.0.36"
        40
    elseif current_version() == v"1.0.37"
        55
    else
        error("Unsupported version: $(current_version())")
    end

    ids_slices = 1:6
    ids_slices_branch = 1:4
    for iter in start_iter:n_sim_iters
        @info "Processing iteration" iter n_sim_iters
        slice_dir = joinpath(data_dir(), "aorta_preprocessed",
                             "v$(current_version().major)", subject, "slices")
        slice_names = "slice_" .* string.(ids_slices)
        slice_branch_names = "slice_branch_" .* string.(ids_slices_branch)
        filenames = vcat(slice_names, slice_branch_names)
        files = joinpath.(slice_dir, filenames) .* ".stl"

        slices = load_geometry.(files)

        results_fluid,
        results_wall = interpolate_fluid_properties(slices, semi; fsi, iter,
                                                    sampling_resolution,
                                                    results_directory,
                                                    wall_thickness,
                                                    interpolate_structure)

        for (i, file) in enumerate(filenames)
            # Plane values
            result_fluid_plane = InitialCondition(;
                                                  particle_spacing=results_fluid[i].sampling_resolution,
                                                  coordinates=results_fluid[i].results_plane.point_coords,
                                                  density=results_fluid[i].results_plane.density,
                                                  pressure=results_fluid[i].results_plane.pressure,
                                                  velocity=results_fluid[i].results_plane.velocity)
            trixi2vtk(result_fluid_plane, filename="fluid_plane_" * file * "_iter_$iter",
                      output_directory=joinpath(results_directory, "slices"),
                      sampling_resolution=sampling_resolution,
                      max_radius=results_fluid[i].max_radius,
                      mean_radius=results_fluid[i].mean_radius,
                      min_radius=results_fluid[i].min_radius)

            result_fluid = InitialCondition(;
                                            particle_spacing=results_fluid[i].sampling_resolution,
                                            coordinates=results_fluid[i].point_coords,
                                            density=results_fluid[i].density,
                                            pressure=results_fluid[i].pressure,
                                            velocity=results_fluid[i].velocity)
            trixi2vtk(result_fluid, filename="fluid_" * file * "_iter_$iter",
                      output_directory=joinpath(results_directory, "slices"),
                      sampling_resolution=sampling_resolution,
                      max_radius=results_fluid[i].max_radius,
                      mean_radius=results_fluid[i].mean_radius,
                      min_radius=results_fluid[i].min_radius)

            # Write center coords
            trixi2vtk(results_fluid[i].center_coords,
                      filename="center_" * file * "_iter_$iter",
                      particle_spacing=0.001 * ones(3),
                      output_directory=joinpath(results_directory, "slices"))

            if fsi && interpolate_structure
                result_wall_plane = InitialCondition(; density=1000.0,
                                                     particle_spacing=results_wall[i].sampling_resolution,
                                                     coordinates=results_wall[i].results_plane.point_coords,
                                                     velocity=results_wall[i].results_plane.velocity)
                trixi2vtk(result_wall_plane, filename="wall_plane_" * file * "_iter_$iter",
                          output_directory=joinpath(results_directory, "slices"),
                          sampling_resolution=sampling_resolution,
                          von_mises_stress=results_wall[i].results_plane.von_mises_stress,
                          displacement=results_wall[i].results_plane.displacement,
                          max_radius=results_wall[i].max_radius,
                          mean_radius=results_wall[i].mean_radius,
                          min_radius=results_wall[i].min_radius)

                result_wall = InitialCondition(; density=1000.0,
                                               particle_spacing=results_wall[i].sampling_resolution,
                                               coordinates=results_wall[i].point_coords,
                                               velocity=results_wall[i].velocity)
                trixi2vtk(result_wall, filename="wall_" * file * "_iter_$iter",
                          output_directory=joinpath(results_directory, "slices"),
                          sampling_resolution=sampling_resolution,
                          von_mises_stress=results_wall[i].von_mises_stress,
                          displacement=results_wall[i].displacement,
                          max_radius=results_wall[i].max_radius,
                          mean_radius=results_wall[i].mean_radius,
                          min_radius=results_wall[i].min_radius)
            end
        end
    end

    @info "VTK slice writing completed" results_directory
    return results_directory
end
