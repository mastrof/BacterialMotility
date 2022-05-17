#==
Periodic boundary conditions have to be passed to a BacterialSystem,
taking only a single AbstractBacterium as argument.
Extra arguments are here provided to allow user customization.
==#

export
    periodic_ortho!,
    reflecting_ortho!,
    sticky_ortho!

function periodic_ortho!(b::AbstractBacterium, xmin::Real, xmax::Real)
    Δx = xmax - xmin
    x = b.r
    for i in eachindex(x)
        if x[i] < xmin
            x[i] += Δx
        elseif x[i] > xmax
            x[i] -= Δx
        end # if
    end # for
end # function

function periodic_ortho!(b::AbstractBacterium, xmin::AbstractVector, xmax::AbstractVector)
    x = b.r
    for i in eachindex(x)
        Δx = xmax[i] - xmin[i]
        if x[i] < xmin[i]
            x[i] += Δx
        elseif x[i] > xmax[i]
            x[i] -= Δx
        end # if
    end # for
end # function

function reflecting_ortho!(b::AbstractBacterium, xmin::Real, xmax::Real)
    x = b.r
    v = b.v
    for i in eachindex(x)
        if x[i] < xmin
            Δx = xmin - x[i]
            x[i] = xmin + Δx
            v[i] *= -1
        elseif x[i] > xmax
            Δx = x[i] - xmax
            x[i] = xmax - Δx
            v[i] *= -1
        end # if
    end # for
end # function

function reflecting_ortho!(b::AbstractBacterium, xmin::AbstractVector, xmax::AbstractVector)
    x = b.r
    v = b.v
    for i in eachindex(x)
        if x[i] < xmin[i]
            Δx = xmin[i] - x[i]
            x[i] = xmin[i] + Δx
            v[i] *= -1
        elseif x[i] > xmax[i]
            Δx = x[i] - xmax[i]
            x[i] = xmax[i] - Δx
            v[i] *= -1
        end # if
    end # for
end # function

function sticky_ortho!(b::AbstractBacterium, xmin::Real, xmax::Real)
    x = b.r
    for i in eachindex(x)
        if x[i] < xmin
            x[i] = xmin
        elseif x[i] > xmax
            x[i] = xmax
        end # if
    end # for
end # function

function sticky_ortho!(b::AbstractBacterium, xmin::AbstractVector, xmax::AbstractVector)
    x = b.r
    for i in eachindex(x)
        if x[i] < xmin[i]
            x[i] = xmin[i]
        elseif x[i] > xmax[i]
            x[i] = xmax[i]
        end # if
    end # for
end # function
