export
    AbstractField, EmptyField,
    AbstractAnalyticField, AnalyticField,
    Concentration_SteadyDiffusionSphericalSource_3D

abstract type AbstractField end

struct EmptyField <: AbstractField
    ;
end # struct

abstract type AbstractAnalyticField <: AbstractField end

@with_kw struct AnalyticField <: AbstractAnalyticField
    ϕ::Function
    ∇ϕ::Function
end # struct

@with_kw struct Concentration_SteadyDiffusionSphericalSource_3D <: AbstractAnalyticField
    R::Float64 # μm, source radius
    C::Float64 # μM, concentration at source surface
    Cbg::Float64 = 0.0 # μM, background concentration
    ϕ::Function = b -> Cbg + C*R/norm(b.r) # concentration field
    ∇ϕ::Function = b -> -C*R/(norm(b.r)^2) # concentration gradient (radial towards center)
end # struct

