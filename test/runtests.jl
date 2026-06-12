using DecoratedParticles
using Test

@testset "DecoratedParticles.jl" begin
    include("test_states.jl")
    include("test_atomsbase.jl")
    include("test_differentiation.jl")
end
