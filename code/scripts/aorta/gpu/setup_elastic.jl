using SimulationSetup
using TrixiParticles
using AMDGPU

parallelization_backend = ROCBackend()
# If `MetalBackend`
# coord_eltype = Float32
coord_eltype = Float64

trixi_include_changeprecision(Float32, @__MODULE__,
                              joinpath(@__DIR__, "..", "setup_elastic.jl"),
                              parallelization_backend=parallelization_backend,
                              coord_eltype=coord_eltype)
