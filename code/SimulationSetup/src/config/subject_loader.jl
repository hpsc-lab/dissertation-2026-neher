"""
Subject-specific parameter loading and extraction.

Provides utilities to load physiological reference data from subject parameter
CSV files and extract commonly used hemodynamic parameters.
"""

"""Load subject-specific parameter row from the parameter CSV.

Retrieves a complete parameter row for a given subject from the external
parameter CSV file, used for subject-specific hemodynamic and material properties.

# Arguments
- `subject::String`: Pseudonym used to look up the row in the CSV.
- `file`: Optional path to the parameter CSV. Defaults to:
          `data_dir()/aorta_initial_condition/parameter_subjects.csv`

Returns:
    DataFrame row filtered by Pseudonym column matching the subject identifier.
"""
function subject_parameters(subject::String;
                            file=joinpath(data_dir(), "aorta_initial_condition",
                                          "parameter_subjects.csv"))::DataFrames.DataFrame
    data = get_parameter_table(; file)
    return data[data.Pseudonym .== subject, :]
end

"""Read the parameter table CSV and return a DataFrame.

Loads the complete subject parameter table from CSV file.

# Arguments
- `file`: Optional path to the parameter CSV. Defaults to:
          `data_dir()/aorta_initial_condition/parameter_subjects.csv`

Returns:
    DataFrame containing all subject parameter records.
"""
function get_parameter_table(;
                             file=joinpath(data_dir(), "aorta_initial_condition",
                                           "parameter_subjects.csv"))::DataFrames.DataFrame
    return read(file, DataFrame)
end

"""Extract commonly used subject parameters as a named tuple.

Extracts key hemodynamic parameters from a subject's parameter row,
applying appropriate unit conversions to SI units.

# Arguments
- `subject::String`: Subject identifier.
- `data`: Optional pre-loaded subject parameters DataFrame row.
          If not provided, fetches via `subject_parameters(subject)`.

Returns:
    Named tuple with fields:
    - `PWV::Float64` - Pulse Wave Velocity (m/s)
    - `A_STJ::Float64` - Sinotubular Junction cross-sectional area (m²)
    - `stroke_volume::Float64` - Cardiac stroke volume (m³)
    - `v_STJ_peak::Float64` - Peak sinotubular junction velocity (m/s)
"""
function extract_subject_parameters(subject::String;
                                    data=subject_parameters(subject))::NamedTuple
    # Pulse Wave Velocity (already in m/s)
    PWV = data[1, Symbol("cfPW Velocity (m/s) carotid-femoral Pulse wave velocity")]

    # Sinotubular Junction Area (convert from cm² to m²)
    A_STJ = data[1, Symbol("STJ_Area (cm^2)")] * cm2_to_m2()

    # Calculate mean blood flow (stroke volume, convert from ml to m³)
    stroke_volume = data[1, Symbol("Schlagvolumen (ml)")] * ml_to_m3()

    # Calculate peak inflow velocity (convert from cm/s to m/s)
    v_STJ_peak = data[1, Symbol("Vmax_STJ (cm/s)")] * cm_to_m()

    return (; PWV, A_STJ, stroke_volume, v_STJ_peak)
end
