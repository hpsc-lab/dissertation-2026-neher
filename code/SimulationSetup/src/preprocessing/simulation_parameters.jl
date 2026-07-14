"""
Simulation parameter assembly and boundary setup.

Constructs SimulationParameters struct by integrating geometric, hemodynamic,
and numerical parameters for a complete cardiovascular simulation setup.
Includes boundary geometry preprocessing and RCR Windkessel model distribution.
"""

"""
    SimulationParameters(subject; particle_spacing, k_d=0.9, p_syst=125.0,
                         p_diast=75.0, L_eff=0.35, T=0.75, density_blood=1060.0,
                         transition_length=transition_length(),
                         boundary_names=boundary_names(),
                         geometry_files=:default)

Assemble subject-specific geometric, hemodynamic, and numerical parameters
for a cardiovascular simulation setup.

This function loads and preprocesses boundary geometries, computes global
lumped-parameter (RCR Windkessel) properties, distributes outlet-specific
Windkessel parameters based on geometry-corrected flow splits, and constructs
all boundary condition data structures required for the simulation.

# Arguments
- `subject::String`: Identifier of the patient used to load geometry
                     and patient-specific reference data.

# Keywords
- `particle_spacing`: Particle spacing used for the simulation.
- `k_d=0.9`: Dimensionless distal resistance fraction used to split the total
             resistance into proximal and distal components.
- `p_syst=125.0`: Target systolic blood pressure in mmHg.
- `p_diast=75.0`: Target diastolic blood pressure in mmHg.
- `L_eff=0.35`: Effective characteristic vessel length used for compliance
                estimation.
- `T=0.75`: Cardiac cycle period.
- `density_blood=1060.0`: Blood density.
- `transition_length`: Axial shift applied to planar boundary geometries.
- `boundary_names`: List of boundary identifiers (including `"inflow"`).
- `geometry_files`: Paths to the boundary geometry files. If `:default`, files are
                    constructed based on the `subject` and `boundary_names`.
"""
struct SimulationParameters{ELTYPE}
    boundaries         :: Dict
    lumped_parametes   :: LumpedParameters
    subject_parameters :: NamedTuple
    T                  :: ELTYPE
    density_blood      :: ELTYPE
    particle_spacing   :: ELTYPE
    n_outlets          :: Int
    subject            :: String

    function SimulationParameters(subject; particle_spacing, k_d=0.9, p_syst=125.0,
                                  p_diast=75.0, L_eff=get_L_eff(subject), T=0.75,
                                  density_blood=1060.0, q_prescribed=realistic_flow_ratios,
                                  boundary_shift=transition_length(),
                                  boundary_names=boundary_names(),
                                  stroke_volume_factor=1.0, v_peak_factor=1.0,
                                  geometry_files=:default)
        ELTYPE = typeof(particle_spacing) # This is for GPUs wich only support Float32

        geometry_files = if geometry_files === :default
            joinpath.(data_dir(), "aorta_preprocessed", "v$(current_version().major)",
                      subject, subject * "_" .* boundary_names .* ".stl")
        else
            geometry_files
        end

        check_outlet_configuration!(geometry_files, boundary_names)

        parameters = extract_subject_parameters(subject)
        (; stroke_volume, A_STJ, v_STJ_peak, PWV) = parameters
        Q_mean = stroke_volume_factor * stroke_volume / T

        # Construct a dictionary mapping boundary names to their corresponding meshes
        geometries = load_geometry.(geometry_files; element_type=ELTYPE)
        geometries_shifted = shift_planar_geometry.(geometries, boundary_shift)
        bnd_geometries = Dict(zip(boundary_names, geometries_shifted))

        rcr_parameters = LumpedParameters(bnd_geometries, Q_mean, density_blood, PWV, A_STJ;
                                          q_prescribed=q_prescribed, k_d=k_d, p_syst=p_syst,
                                          p_diast=p_diast, L_eff=L_eff)
        (; R_T, C_T, k_d, q) = rcr_parameters
        pressure_model_dict = Dict{String, RCRWindkesselModel}()
        check_sum_q = 0.0
        for key in keys(bnd_geometries)
            key == "inflow" && continue
            R_1 = (1 - k_d) * R_T / q[key] # R_1,i = 0.1 * p_mean / Q_i = p_mean / (q_i * Q_mean)
            R_2 = k_d * R_T / q[key] # R_2,i = 0.9 * p_mean / Q_i = p_mean / (q_i * Q_mean)
            C = C_T * q[key] # C,i = C_T * q_i (C = ∑C_i)
            check_sum_q += q[key]
            pressure_model_dict[key] = RCRWindkesselModel(; characteristic_resistance=R_1,
                                                          peripheral_resistance=R_2,
                                                          compliance=C)
        end

        if !isapprox(check_sum_q, 1.0)
            @warn "Outlet flow fractions sum to $check_sum_q (expected 1.0)"
        end

        inf_ = convert(ELTYPE, Inf) # This is for GPUs wich only support Float32
        boundary_dict = Dict("inflow" => (geometry=bnd_geometries["inflow"], id=0,
                                          transition_length=boundary_shift,
                                          sample_points=create_sample_points(bnd_geometries["inflow"],
                                                                             particle_spacing),
                                          pressure_model=RCRWindkesselModel(;
                                                                            characteristic_resistance=inf_,
                                                                            peripheral_resistance=inf_,
                                                                            compliance=inf_),
                                          cross_sectional_area=surface_area(bnd_geometries["inflow"])))
        outlet_ids = boundary_ids()
        for key in keys(bnd_geometries)
            key == "inflow" && continue
            boundary_dict[key] = (geometry=bnd_geometries[key], id=outlet_ids[key],
                                  transition_length=boundary_shift,
                                  sample_points=create_sample_points(bnd_geometries[key],
                                                                     particle_spacing),
                                  pressure_model=pressure_model_dict[key],
                                  cross_sectional_area=surface_area(bnd_geometries[key]))
        end

        n_outlets = haskey(boundary_dict, "brachiocephalic") ? 4 : 5

        # Calculate peak inflow velocity
        A_in = surface_area(boundary_dict["inflow"].geometry)

        v_peak = v_peak_factor * v_STJ_peak * (A_STJ / A_in)

        data_subject = subject_parameters(subject)

        parameters = (; Q_mean, A_in, v_peak=convert(ELTYPE, v_peak),
                      E=convert(ELTYPE, data_subject[1, "Young's Modulus (Mpa)"] * 1e6),
                      nu=convert(ELTYPE, data_subject[1, "Poisson's ratio (v)"]),
                      p_syst=p_syst, p_diast=p_diast, L_eff=L_eff,
                      p_mean=(p_syst + 2 * p_diast) / 3, data=data_subject, parameters...)

        return new{ELTYPE}(boundary_dict, rcr_parameters, parameters, T, density_blood,
                           particle_spacing, n_outlets, subject)
    end
end

@inline Base.eltype(::SimulationParameters{ELTYPE}) where {ELTYPE} = ELTYPE

function Base.show(io::IO, ::MIME"text/plain", obj::SimulationParameters)
    @nospecialize obj # reduce precompilation time

    if get(io, :compact, false)
        show(io, obj)
    else
        summary_header(io, "SimulationParameters{$(eltype(obj))}")
        summary_line(io, "subject", obj.subject)
        summary_line(io, "#outlets", obj.n_outlets)
        summary_line(io, "density_blood", "$(obj.density_blood) kg/m³")
        summary_line(io, "T", "$(obj.T) sec. ($(60/obj.T) bpm)")
        summary_line(io, "Q_mean", "$(obj.subject_parameters.Q_mean * m3_to_ml()) ml/sec")
        summary_line(io, "Mean Pressure", "$(obj.subject_parameters.p_mean) mmHg")
        summary_line(io, "p_diast", "$(obj.subject_parameters.p_diast) mmHg")
        summary_line(io, "p_syst", "$(obj.subject_parameters.p_syst) mmHg")
        summary_line(io, "Velocities", "")
        summary_line(io, "  Inflow Peak", "$(obj.subject_parameters.v_peak) m/s")
        summary_line(io, "  STJ Peak", "$(obj.subject_parameters.v_STJ_peak) m/s")
        summary_line(io, "Areas", "")
        summary_line(io, "  Inflow", "$(obj.subject_parameters.A_in * m2_to_cm2()) cm²")
        summary_line(io, "  STJ", "$(obj.subject_parameters.A_STJ * m2_to_cm2()) cm²")
        summary_line(io, "E-Modulus", "$(obj.subject_parameters.E * 1e-6) MPa")
        summary_line(io, "Poisson's ratio", "$(obj.subject_parameters.nu)")
        summary_footer(io)
    end
end

function pad_planar_face(face, padding)
    padding == 0 && return face

    v1, v2, v3 = (face[:, 1], face[:, 2], face[:, 3])
    edge1 = v2 - v1
    edge2 = v3 - v1
    dir1 = normalize(edge1)
    dir2 = normalize(edge2)

    v1_padded = v1 - padding * dir1 - padding * dir2
    v2_padded = v2 + padding * dir1 - padding * dir2
    v3_padded = v3 - padding * dir1 + padding * dir2

    return stack((v1_padded, v2_padded, v3_padded))
end

function create_sample_points(geometry::TriangleMesh, particle_spacing; clip=true,
                              padding=zero(particle_spacing))
    boundary_fac_, face_normal = planar_geometry_to_face(geometry)
    if padding > 0
        face = pad_planar_face(stack(boundary_fac_), padding)
        boundary_fac = (face[:, 1], face[:, 2], face[:, 3])
    else
        boundary_fac = boundary_fac_
    end

    ic = @suppress extrude_geometry(boundary_fac; particle_spacing, density=1,
                                    direction=face_normal, n_extrude=1)
    !(clip) && return ic.coordinates

    cutter = extrude_geometry(geometry, particle_spacing * 3 / 2)

    return intersect(ic, cutter).coordinates
end

function create_sample_points(point_cloud::AbstractArray, resolution;
                              padding=zero(resolution),
                              intersection_mask=nothing, particle_spacing)
    face_vertices_ = TrixiParticles.oriented_bounding_box(point_cloud)
    face_vertices = pad_planar_face(face_vertices_, padding)

    # Vectors spanning the face
    edge1 = face_vertices[:, 2] - face_vertices[:, 1]
    edge2 = face_vertices[:, 3] - face_vertices[:, 1]

    face_normal = SVector(Tuple(normalize(cross(edge1, edge2))))

    face = (face_vertices[:, 1], face_vertices[:, 2], face_vertices[:, 3])

    ic = @suppress extrude_geometry(face; particle_spacing=resolution, density=1,
                                    direction=face_normal, n_extrude=1)

    isnothing(intersection_mask) && return ic.coordinates

    too_close = find_too_close_particles(ic.coordinates, intersection_mask,
                                         particle_spacing)

    return ic.coordinates[:, too_close]
end

function check_outlet_configuration!(files, geometry_boundaries)
    missing_idx = findall(f -> !isfile(f), files)
    if Set(geometry_boundaries[missing_idx]) == Set(["right_subclavian", "right_common"]) ||
       Set(geometry_boundaries[missing_idx]) == Set(["brachiocephalic"])
        deleteat!(files, missing_idx)
        deleteat!(geometry_boundaries, missing_idx)
    elseif isempty(missing_idx)
        error("To much outlet files")
    else
        error("Unexpected missing outlet files: $(geometry_boundaries[missing_idx])")
    end
end

function check_pressure_models!(open_boundary, boundary_dict)
    if open_boundary.cache.pressure_reference_values[2] ==
       boundary_dict["thoracic"].pressure_model &&
       open_boundary.cache.pressure_reference_values[3] ==
       boundary_dict["left_common"].pressure_model &&
       open_boundary.cache.pressure_reference_values[4] ==
       boundary_dict["left_subclavian"].pressure_model
    else
        error("boundary zones are not sorted")
    end

    if haskey(boundary_dict, "brachiocephalic") &&
       open_boundary.cache.pressure_reference_values[5] ==
       boundary_dict["brachiocephalic"].pressure_model
    elseif haskey(boundary_dict, "right_subclavian") &&
           haskey(boundary_dict, "right_common") &&
           open_boundary.cache.pressure_reference_values[5] ==
           boundary_dict["right_subclavian"].pressure_model &&
           open_boundary.cache.pressure_reference_values[6] ==
           boundary_dict["right_common"].pressure_model
    else
        error("unkown configuration")
    end
end

function get_L_eff(subject)
    file_aorta = joinpath(data_dir(), "aorta_preprocessed", "v$(current_version().major)",
                          subject, subject * ".stl")
    geometry_aorte = load_geometry(file_aorta)

    return norm(geometry_aorte.min_corner - geometry_aorte.max_corner)
end
