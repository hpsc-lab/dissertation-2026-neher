using SimulationSetup
using TrixiParticles
using Metal

parallelization_backend = MetalBackend()
coord_eltype = parallelization_backend isa MetalBackend ? Float32 : Float64

trixi_include_changeprecision(Float32, @__MODULE__,
                              joinpath(@__DIR__, "..", "setup_rigid.jl"),
                              parallelization_backend=parallelization_backend,
                              coord_eltype=coord_eltype)
