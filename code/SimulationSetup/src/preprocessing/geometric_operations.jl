"""Shift a planar geometry along its averaged face normal.

Returns a TriangleMesh translated by `shift_length` in the direction of the
average face normal. Useful to create offset boundary surfaces.
"""
function shift_planar_geometry(geometry, shift_length)
    face_normal = normalize(sum(geometry.face_normals) / nfaces(geometry))
    shift = face_normal * shift_length

    vertices_origin = copy(geometry.vertices)
    faces_origin = copy(geometry.face_vertices)
    normals_shifted = copy(geometry.face_normals)

    # Shift every origin vertex by the same vector to create the shifted vertices.
    vertices_shifted = [v .+ shift for v in vertices_origin]
    # Shift each face-tuple component-wise to create shifted face tuples
    face_vertices_shifted = [(v1 .+ shift, v2 .+ shift, v3 .+ shift)
                             for (v1, v2, v3) in faces_origin]

    return TriangleMesh(face_vertices_shifted, normals_shifted, vertices_shifted)
end

"""Compute the surface area of a planar triangle mesh.

Uses the cross-product formula for triangle area and sums over faces.
"""
function surface_area(mesh::TriangleMesh)
    area = sum(mesh.face_vertices) do vertices
        # Formula for the area of a triangle using cross product:
        # A = ||(v2 - v1) × (v3 - v1)|| / 2
        v1, v2, v3 = vertices
        edge1 = v2 - v1
        edge2 = v3 - v1
        return norm(cross(edge1, edge2)) / 2
    end

    return area
end

"""Project points orthogonally onto a plane defined by `rectangular_face` and its normal.

Arguments:
- points: 3×N array of point coordinates
- face_normal: normal vector of the plane
- rectangular_face: collection of corner points; first corner is used as origin
"""
function project_points_to_plane(points, face_normal, rectangular_face)
    n = normalize(face_normal)

    # Plane origin (first corner point of the rectangle)
    p0 = rectangular_face[1]

    # Orthogonal projection
    d = points .- p0  # Offset from p0, 3×N
    dist_to_plane = vec(n' * d)  # Distance in normal direction, N

    # Projected points: points - dist_to_plane[i] * n for each point i
    projected = points .- n .* reshape(dist_to_plane, 1, :)

    return projected
end
