export AbstractBacterium, AbstractBacterialSystem

@doc raw"""
    AbstractBacterium{D,T}

General interface for bacterium object.
    - `D`: dimensionality of the space in which the bacterium lives
    - `T`: type of the space in which the bacterium lives
"""
abstract type AbstractBacterium{D,T} end


abstract type AbstractBacterialSystem end

