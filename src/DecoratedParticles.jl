module DecoratedParticles

# core of the package : manipulating particle states
include("states.jl")

# ForwardDiff differentiation w.r.t. states and NamedTuples
include("differentiation.jl")

include("show.jl")

include("atomsbase.jl")


end
