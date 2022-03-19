export
    run!,
    rotate2D!, rotate3D!,
    Degenerate,
    tumble!, reverse!, flick!, reverse_flick!


@doc raw"""
    run!(bacterium::B, dt) where B<:AbstractBacterium

Update bacterial position after a ballistic run of length dt.
"""
function run!(bacterium::B, dt) where B<:AbstractBacterium
    @. bacterium.r += bacterium.v * dt
end # function

@doc raw"""
    rotate2D!(w, θ)

Rotate two-dimensional vector w of an angle θ.
"""
function rotate2D!(w, θ)
    w .= Angle2d(θ) * w
end # function

@doc raw"""
    rotate3D!(w, θ)

Rotate three-dimensional vector w of an angle θ (in radians) with uniform off-axis distribution.
"""
function rotate3D!(w, θ)
    m = findfirst(w .≠ 0.0)
    n = m%3 + 1
    u = zero(w)
    u[n] = w[m]
    u[m] = -w[n]
    a = AngleAxis(θ, u...) * w
    ϕ = rand() * 2π
    w .= AngleAxis(ϕ, w...) * a
end # function


@doc raw"""
    tumble!(bacterium::B, PDF=Uniform(0,2π)) where B<:AbstractBacterium{D}

Bacterium performs a tumble in D-dimensional space with angular distribution given by PDF.
Implemented for D=2 and D=3.
"""
tumble!(bacterium::B, PDF=Uniform(0,2π)) where B<:AbstractBacterium{2} = rotate2D!(bacterium.v, rand(PDF))
tumble!(bacterium::B, PDF=Uniform(0,2π)) where B<:AbstractBacterium{3} = rotate3D!(bacterium.v, rand(PDF))


struct Degenerate{T<:Real} <: ContinuousUnivariateDistribution
    x::T
end # struct

Base.rand(d::Degenerate) = d.x

reverse!(bacterium::B, PDF=Degenerate(π)) where B<:AbstractBacterium{2} = rotate2D!(bacterium.v, rand(PDF))
reverse!(bacterium::B, PDF=Degenerate(π)) where B<:AbstractBacterium{3} = rotate3D!(bacterium.v, rand(PDF))
flick!(bacterium::B, PDF=Degenerate(π/2)) where B<:AbstractBacterium{2} = rotate2D!(bacterium.v, rand(PDF))
flick!(bacterium::B, PDF=Degenerate(π/2)) where B<:AbstractBacterium{3} = rotate3D!(bacterium.v, rand(PDF))

@doc raw"""
    reverse_flick!(bacterium::B, PDFreverse=Degenerate(π), PDFflick=Degenerate(π/2)) where B<:AbstractBacterium{D}

Depending on its internal `RunState` bacterium performs either a reverse or a flick, with angular distributions respectively PDFreverse and PDFflick.
If the bacterium does not have a `RunState` entry as state variable, it is created in place and initialized randomly to either 1 (reverse) or 0 (flick).
After a reverse, the run state is switched to flick mode and viceversa.
Implemented for D=2 and D=3.
"""
function reverse_flick!(bacterium::B, PDFreverse=Degenerate(π), PDFflick=Degenerate(π/2)) where B<:AbstractBacterium
    # if bacterium doesn't have a "RunState" key, create it with random initial value 0 or 1
    if !haskey(bacterium.state, "RunState")
        bacterium.state["RunState"] = rand([0,1])
    end # if
    if Bool(bacterium.state["RunState"]) # == 1 then reverse
        reverse!(bacterium, PDFreverse)
        bacterium.state["RunState"] = 0
    else # == 0 then flick
        flick!(bacterium, PDFflick)
        bacterium.state["RunState"] = 1
    end # if
end # function
