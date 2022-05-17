module BacterialMotility

using StaticArrays
using LinearAlgebra
using Distributions
using Rotations
using Parameters

# abstract types
include("interface.jl")

# utility functions
include("utils.jl")

# core
include("properties.jl")
include("fields.jl")
include("motion.jl")
include("sensing.jl")
include("bacterium.jl")
include("bacterialsystem.jl")

# implementations of chemotaxis models
include("Brumley.jl")
include("BrownBerg.jl")

end # module
