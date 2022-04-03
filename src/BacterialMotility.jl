module BacterialMotility

using StaticArrays
using LinearAlgebra
using Distributions
using Rotations
using Parameters


function dummy(args...; kwargs...)
    ;
end

struct Degenerate{T<:Real} <: ContinuousUnivariateDistribution
    x::T
end # struct

Base.rand(d::Degenerate) = d.x


include("properties.jl")
include("bacterium.jl")
include("fields.jl")
include("motion.jl")
include("sensing.jl")
include("integration.jl")

# model-specific implementations
include("Brumley.jl")
include("BrownBerg.jl")

end # module
