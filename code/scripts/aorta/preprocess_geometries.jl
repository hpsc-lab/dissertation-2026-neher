# Preprocessing script for aorta geometries.

# Loads STL geometries from the segmentation dataset, applies preprocessing
# steps (centering and scaling), and saves transformed geometries.

# Steps:
#     1. Load STL geometries for cases F01-F16 (excluding F04 and F06)
#     2. Center geometries at the origin
#     3. Scale by factor 1e-3 (mm to m conversion)
#     4. Save preprocessed geometries to aorta_centered directory
using SimulationSetup

is_in_dataset = setdiff(1:16, (4, 6))
geometry_names = ["F$(lpad(i,2,'0')).stl" for i in is_in_dataset]

output_directory = joinpath(data_dir(), "aorta_centered")
isdir(output_directory) || mkdir(output_directory)

for geometry_name in geometry_names
    file = joinpath(data_dir(), "segmentation_aorta_UKA",
                    "Segmentation_Aorta_" * geometry_name)

    geometry = SimulationSetup.load_geometry(file)

    center = (geometry.min_corner .+ geometry.max_corner) ./ 2

    # Center geometry
    map!(v -> v - center, geometry.vertices, geometry.vertices)
    map!(fv -> (fv[1] - center, fv[2] - center, fv[3] - center),
         geometry.face_vertices, geometry.face_vertices)

    # Scale geometry
    map!(v -> v .* 1e-3, geometry.vertices, geometry.vertices)
    map!(fv -> (fv[1] .* 1e-3, fv[2] .* 1e-3, fv[3] .* 1e-3),
         geometry.face_vertices, geometry.face_vertices)

    save_stl(joinpath(output_directory, geometry_name), geometry)
end
