"""
Hemodynamic parameters and lumped-parameter models.

Defines the LumpedParameters struct for RCR Windkessel model properties and
provides various flow ratio strategies for outlet boundary condition setup.
"""

"""
    LumpedParameters(geometries, Q_mean, density_blood, PWV, A_STJ;
                     k_d=0.9, p_syst=125.0, p_diast=75.0, L_eff=0.35)

Compute global lumped-parameter for the `RCRWindkesselModel` properties from geometric
and hemodynamic reference data.
It derives geometry-corrected outlet flow fractions by transferring reference
physiological flow splits to the current geometry using outlet surface areas.

# Arguments
- `geometries`: Dictionary mapping boundary names to boundary geometries (`TriangleMesh`),
                representing the inflow and all outflow surfaces. The key `"inflow"` is
                ignored when computing outlet flow splits.
- `Q_mean`: Mean volumetric inflow rate.
- `density_blood`:Blood density.
- `PWV`: Pulse wave velocity, used to estimate total arterial compliance
          via the Bramwell-Hill relation.
- `A_STJ`: Cross-sectional area at the sinotubular junction (STJ), used as a reference
           area for compliance estimation.

# Keywords
- `k_d=0.9`: Dimensionless distal resistance fraction, used to split the total resistance into
             proximal and distal components. Default: `0.9`.
- `p_syst=125.0`: Target systolic blood pressure in mmHg. Default: `125.0`.
- `p_diast=75.0`: Target diastolic blood pressure in mmHg. Default: `75.0`.
- `L_eff=0.35`: Effective characteristic vessel length used in the compliance
                estimation. Default: `0.35`.
"""
struct LumpedParameters{ELTYPE}
    R_T :: ELTYPE      # Total resistance
    C_T :: ELTYPE      # Total compliance
    k_d :: ELTYPE      # Distal resistance fraction
    q   :: Dict{String, ELTYPE}  # Outlet flow fractions

    function LumpedParameters(geometries, Q_mean, density_blood, PWV, A_STJ;
                              q_prescribed=realistic_flow_ratios,
                              k_d=0.9, p_syst=125.0, p_diast=75.0, L_eff=0.35)
        ELTYPE = typeof(density_blood)

        p_syst_ = p_syst * mmHg_to_Pa()
        p_diast_ = p_diast * mmHg_to_Pa()

        p_mean = (p_syst_ + 2 * p_diast_) / 3
        R_T = convert(ELTYPE, p_mean / Q_mean)

        # Bramwell-Hill-Beziehung
        C_T = convert(ELTYPE, L_eff * A_STJ / (density_blood * PWV^2))

        return new{typeof(R_T)}(R_T, C_T, k_d, q_prescribed(geometries))
    end
end

@inline Base.eltype(::LumpedParameters{ELTYPE}) where {ELTYPE} = ELTYPE

function Base.show(io::IO, ::MIME"text/plain", obj::LumpedParameters)
    @nospecialize obj # reduce precompilation time

    if get(io, :compact, false)
        show(io, obj)
    else
        summary_header(io, "LumpedParameters{$(eltype(obj))}")
        summary_line(io, "Total Resistance (R_T)", "$(obj.R_T) Pa·s/m³")
        summary_line(io, "Total Compliance (C_T)", "$(obj.C_T) m³/Pa")
        summary_line(io, "Distal Resistance Fraction", obj.k_d)
        summary_line(io, "Outlet Flow Fractions (q)", "")
        for (key, value) in obj.q
            summary_line(io, "  $key", "$(round(value*100, digits=3)) %")
        end
        summary_footer(io)
    end
end

"""Compute geometry-corrected flow fractions using reference distributions.

The method transfers physiologic reference flow splits to the current geometry
by scaling with outlet surface areas and normalizing the result.
"""
function baeumler_flow_ratios(geometries)
    ELTYPE = eltype(geometries["inflow"])
    # Note:
    # Bäumler et al., 2024 (supplementary material):
    # "Three-element windkessel parameters were tuned to match the patient's blood pressure
    # of 125/75 mmhg for the individual models."
    #
    # ------------------------------------------------------------------------------
    # FLOW-SPLIT ADJUSTMENT USING REFERENCE VALUES + GEOMETRIC CORRECTION
    #
    # Step 1 — Compute reference outlet velocities
    # -------------------------------------------
    # The reference flow ratios q_ref represent physiologically expected relative
    # flow distributions in a "reference geometry" with outlet areas A_ref_out.
    # In that reference geometry, the mean velocities satisfy:
    #
    #     v_ref[i] ∝ q_ref[i] / A_ref_out[i]
    #
    # (Absolute flow rate cancels out — only relative patterns matter.)
    #
    #
    # Step 2 — Transfer these reference velocities to the new geometry
    # ----------------------------------------------------------------
    # We assume that the characteristic velocity pattern across outlets remains
    # approximately the same in the new patient geometry. Therefore:
    #
    #     v_new[i] ∝ v_ref[i]  ∝ q_ref[i] / A_ref_out[i]
    #
    # Using the new outlet areas A_out[i], the corresponding flow rates become:
    #
    #     Q_new[i] = v_new[i] * A_out[i]
    #              ∝ q_ref[i] * (A_out[i] / A_ref_out[i])
    #
    #
    # Step 3 — Normalize to obtain updated flow fractions
    # ---------------------------------------------------
    # The new, geometry-corrected flow distribution q_new is obtained by
    # normalizing the Q_new[i] so that the sum over all outlets is 1:
    #
    #     q_new[i] =
    #         ( q_ref[i] * (A_out[i] / A_ref_out[i]) ) /
    #         Σ_j ( q_ref[j] * (A_out[j] / A_ref_out[j]) )
    #
    # This preserves the physiological meaning of q_ref while accounting for
    # geometric deviations between the reference model and the current geometry.

    # Calculate the target flow ratio per outlet (q_total = ∑q_i = 1.0).
    # Volume flow ratios reported in the paper:
    # q_BCT = 0.15, q_LCA = 0.031, q_LSA = 0.054, q_outlet = 0.765
    #
    # files = joinpath.(DATA_DIR, "aorta_dissection_stanford", "reference_outlets",
    #                   "aorta_ref_" .* geometry_boundaries .* ".stl")
    # geometries_ref = load_geometry.(files)
    # A_ref_out = surface_area.(geometries_ref)
    # diameters_out_ref = sqrt.(A_ref_out ./ π) .* 2 .* 1000 # mm
    q_ref = Dict(
        "right_subclavian" => (q=0.075, A=6.10542213559623e-5), # A_brachiocephalic / 2
        "right_common" => (q=0.075, A=6.10542213559623e-5), # A_brachiocephalic / 2
        "brachiocephalic" => (q=0.15, A=0.0001221084427119246),
        "left_common" => (q=0.031, A=3.82946155006755e-5),
        "left_subclavian" => (q=0.054, A=5.002637187075045e-5),
        "thoracic" => (q=0.765, A=0.0003443689946634431)
    )

    weights = Dict{String, Float64}()
    for key in keys(geometries)
        key == "inflow" && continue
        A = surface_area(geometries[key])
        weights[key] = q_ref[key].q * A / q_ref[key].A
    end

    w_sum = sum(values(weights))

    return Dict(key => convert(ELTYPE, w / w_sum) for (key, w) in weights)
end

"""Return prescribed realistic flow fractions for typical patient anatomy.

Returns a Dict mapping outlet names to normalized flow fractions. The returned
values are converted to the element type of the inflow geometry.
"""
function realistic_flow_ratios(geometries)
    ELTYPE = eltype(geometries["inflow"])
    prescribed_ratios = Dict(
        "right_subclavian" => (q = 0.075),
        "right_common" => (q = 0.075),
        "brachiocephalic" => (q = 0.15),
        "left_common" => (q = 0.04),
        "left_subclavian" => (q = 0.06),
        "thoracic" => (q = 0.75)
    )

    return Dict(key => convert(ELTYPE, q) for (key, q) in prescribed_ratios)
end

"""Return an intentionally skewed/unrealistic flow split for testing.

Useful for sensitivity testing where outlets receive non-physiologic flow fractions.
The returned values are converted to the element type of the inflow geometry.
"""
function unrealistic_flow_ratios(geometries)
    ELTYPE = eltype(geometries["inflow"])
    prescribed_ratios = Dict(
        "right_subclavian" => (q = 0.125),
        "right_common" => (q = 0.125),
        "brachiocephalic" => (q = 0.25),
        "left_common" => (q = 0.25),
        "left_subclavian" => (q = 0.25),
        "thoracic" => (q = 0.25)
    )

    return Dict(key => convert(ELTYPE, q) for (key, q) in prescribed_ratios)
end

function adapted_flow_ratios(geometries)
    ELTYPE = eltype(geometries["inflow"])
    prescribed_ratios = Dict(
        "right_subclavian" => (q = 0.08),
        "right_common" => (q = 0.1),
        "brachiocephalic" => (q = 0.18),
        "left_common" => (q = 0.075),
        "left_subclavian" => (q = 0.095),
        "thoracic" => (q = 0.65)
    )

    return Dict(key => convert(ELTYPE, q) for (key, q) in prescribed_ratios)
end
