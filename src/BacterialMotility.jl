module BacterialMotility

using StaticArrays
using LinearAlgebra
using Distributions
using Rotations
using Parameters

include("properties.jl")
include("types.jl")
include("fields.jl")
include("motion.jl")
include("sensing.jl")
include("integration.jl")

# model-specific implementations
include("Brumley.jl")

end # module
