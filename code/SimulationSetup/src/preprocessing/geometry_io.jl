"""Helpers to write TriangleMesh objects to binary STL files.

Functions:
- save_stl(filename, mesh; faces=...): convenience wrapper that accepts a
  filename string.
- save_stl(fn::File, mesh; faces=...): accepts a TrixiParticles.File wrapper.
- save_stl(f::Stream, mesh; faces=...): low-level writer that writes binary
  STL according to the STL spec.
"""

function save_stl(filename, mesh; faces=TrixiParticles.eachface(mesh))
    save_stl(TrixiParticles.FileIO.File{TrixiParticles.FileIO.format"STL_BINARY"}(filename),
             mesh; faces)
end

function save_stl(fn::TrixiParticles.FileIO.File{TrixiParticles.FileIO.format"STL_BINARY"},
                  mesh::TrixiParticles.TriangleMesh; faces=TrixiParticles.eachface(mesh))
    open(fn, "w") do s
        save_stl(s, mesh; faces)
    end
end

function save_stl(f::TrixiParticles.FileIO.Stream{TrixiParticles.FileIO.format"STL_BINARY"},
                  mesh::TrixiParticles.TriangleMesh; faces)
    io = TrixiParticles.FileIO.stream(f)
    points = mesh.face_vertices
    normals = mesh.face_normals

    # Implementation follows the binary STL format (80-byte header, uint32 count, triangles)
    for i in 1:80 # Write empty header
        write(io, 0x00)
    end

    write(io, UInt32(length(faces))) # Write triangle count
    for i in faces
        n = SVector{3, Float32}(normals[i])
        triangle = points[i]

        # Normal vector (3 floats)
        for j in 1:3
            write(io, n[j])
        end

        # Three triangle vertices (each 3 floats)
        for point in triangle, p in point
            write(io, Float32(p))
        end
        write(io, 0x0000) # 16-bit attribute byte count (commonly zero)
    end
end
