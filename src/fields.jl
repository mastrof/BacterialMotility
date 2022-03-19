export
    AbstractField, EmptyField, AnalyticField,
    Concentration_SteadyDiffusionSphericalSource_3D

abstract type AbstractField end

struct EmptyField <: AbstractField
    ;
end # struct

@with_kw struct AnalyticField <: AbstractField
    ϕ::Function
    ∇ϕ::Function
end # struct

@with_kw struct Concentration_SteadyDiffusionSphericalSource_3D <: AbstractField
    R::Float64 # μm, source radius
    C::Float64 # μM, concentration at source surface
    Cbg::Float64 = 0.0 # μM, background concentration
    c = r -> Cbg + C*R/r # concentration field
    ∇c = r -> -C*R/(r*r) # concentration gradient (radial towards center)
end # struct
    
