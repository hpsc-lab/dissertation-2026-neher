"""
Boundary configuration and identifiers.

Defines the standard aorta boundary nomenclature and associated identifiers
used across the codebase for consistent indexing and reference.
"""

"""Return a dictionary mapping boundary name strings to numeric IDs.

Used to index outlet arrays consistently across the codebase.

Returns:
    Dict with entries:
    - "inflow" => 0
    - "thoracic" => 1
    - "left_common" => 2
    - "left_subclavian" => 3
    - "brachiocephalic" => 4
    - "right_subclavian" => 4
    - "right_common" => 5
"""
function boundary_ids()::Dict{String, Int}
    return Dict("inflow" => 0,
                "thoracic" => 1,
                "left_common" => 2,
                "left_subclavian" => 3,
                "brachiocephalic" => 4,
                "right_subclavian" => 4,
                "right_common" => 5)
end

"""Return the ordered list of boundary names used by the package.

The ordering is significant for geometric preprocessing and parameter assignment.

Returns:
    Vector of strings: ["inflow", "thoracic", "left_common", "left_subclavian",
                        "brachiocephalic", "right_subclavian", "right_common"]
"""
function boundary_names()::Vector{String}
    return ["inflow",
        "thoracic",
        "left_common",
        "left_subclavian",
        "brachiocephalic",
        "right_subclavian",
        "right_common"]
end
