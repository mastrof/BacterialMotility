export
    run!, rotate!,
    Degenerate,
    rotational_diffusion!,
    tumble!, reverse!, flick!, reverse_flick!


@doc raw"""
    run!(bacterium::AbstractBacterium, dt)

Update bacterial position after a ballistic run of length dt.
"""
function run!(bacterium::AbstractBacterium, dt)
    @. bacterium.r += bacterium.v * dt
end # function

function rotate!(w::T) where T<:StaticVector{1}
     w .= -w
end # function

function rotate!(w::T, θ) where T<:StaticVector{2}
    w .= Angle2d(θ) * w
end # function

function rotate!(w::T, θ) where T<:StaticVector{3}
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
    tumble!(bacterium::B, PDF=Uniform(0,2π)) where B<:AbstractBacterium

Bacterium performs a tumble in D-dimensional space with angular distribution given by PDF.
Implemented for D=1,2,3. In D=1 the tumble is just a reversal.
"""
tumble!(bacterium::AbstractBacterium{1}, PDF=Uniform(0,2π)) = rotate!(bacterium.v)
tumble!(bacterium::AbstractBacterium, PDF=Uniform(0,2π)) = rotate!(bacterium.v, rand(PDF))


@doc raw"""
    rotational_diffusion!(bacterium::AbstractBacterium, PDF=σ->Normal(0,σ))
Perform random rotation of `bacterium` due to rotational diffusion. The turn angle is chosen from a normal distribution whose standard deviation is given by `σ = sqrt(2*Drot*dt)` and `Drot` and `dt` are specified in `bacterium.state`.
If `bacterium` is an AbstractBacterium{1}, the function does nothing.
"""
function rotational_diffusion!(bacterium::AbstractBacterium,
                               PDF=σ->Normal(0,σ))
    Drot = bacterium.state["RotationalDiffusivity"] # rad²/s
    dt = bacterium.state["IntegrationTimestep"] # s
    σ = sqrt(2*Drot*dt) # rad
    rotate!(bacterium.v, rand(PDF(σ)))
end # function

function rotational_diffusion!(bacterium::AbstractBacterium{1},
                               PDF=σ->Normal(0,σ))
    nothing
end # function


reverse!(bacterium::AbstractBacterium, PDF=Degenerate(π)) = rotate!(bacterium.v, rand(PDF))
flick!(bacterium::AbstractBacterium, PDF=Degenerate(π/2)) = rotate!(bacterium.v, rand(PDF))

@doc raw"""
    reverse_flick!(bacterium::AbstractBacterium, PDFreverse=Degenerate(π), PDFflick=Degenerate(π/2))

Depending on its internal `RunState` bacterium performs either a reverse or a flick, with angular distributions respectively PDFreverse and PDFflick.
After a reverse, the run state is switched to flick mode and viceversa.
Implemented for D=2 and D=3.
"""
function reverse_flick!(bacterium::AbstractBacterium, PDFreverse=Degenerate(π), PDFflick=Degenerate(π/2))
    s = bacterium.state["RunState"]
    if s == 1 # then reverse
        reverse!(bacterium, PDFreverse)
    else # == 0 then flick
        flick!(bacterium, PDFflick)
    end # if
    bacterium.state["RunState"] = 1-s
end # function
