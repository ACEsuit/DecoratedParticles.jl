

using DecoratedParticles, StaticArrays, Test, LinearAlgebra
using DecoratedParticles: PState, VState
DP = DecoratedParticles

##

X = PState( rr = randn(SVector{3, Float64}) )

SYMS = DP._syms(X)
println(@test SYMS == (:rr,))

TT = DP._tt(X)
println(@test( TT == Tuple{SVector{3, Float64}} ))

TDX = DP.vstate_type(X)
println(@test TDX == VState{NamedTuple{(:rr,), Tuple{SVector{3, Float64}}}})


dX = VState(X)
println(@test typeof(dX) == DP.vstate_type(X))
println(@test dX.rr == X.rr)


cTDX = DP.vstate_type(0.0im, X)
println(@test cTDX == VState{NamedTuple{(:rr,), Tuple{SVector{3, ComplexF64}}}})

cdX = complex(dX)
println(@test( cdX == VState(rr = X.rr .+ 0im) ))
println(@test( real(cdX) == dX ))
println(@test( imag(cdX) == VState(rr = zero(SVector{3, Float64})) ))

# not sure how to test this, but at least it should work:
@show rand(TDX)
@show randn(TDX)
@show zero(TDX)

@show rand(X)
@show randn(X)
@show zero(X)

@info("arithmetic")
println(@test( X + dX == PState(rr = X.rr + dX.rr) ))
println(@test( X - dX == PState(rr = X.rr - dX.rr) ))
println(@test( dX + cdX == VState(rr = dX.rr + cdX.rr) ))
println(@test( dX - cdX == VState(rr = dX.rr - cdX.rr) ))

a = randn() 
println(@test( a * dX == VState(rr = a * dX.rr) ))
println(@test( dX * a == VState(rr = a * dX.rr) ))
println(@test( -dX == VState(rr = -dX.rr) ))

println(@test( dot(dX, cdX) == dot(dX.rr, cdX.rr) ))
println(@test( dot(cdX, dX) == dot(cdX.rr, dX.rr) ))

dX1, dX2 = randn(cdX), randn(cdX)
DP.contract(dX1, dX2)
println(@test DP.contract(dX1, dX2) == sum(dX1.rr .* dX2.rr) )
println(@test DP.contract(dX1, dX2) != dot(dX1, dX2) )

println(@test isapprox(dX, cdX))

println(@test norm(dX) == norm(dX.rr) )
println(@test DP.sumsq(dX) == DP.sumsq(dX.rr) )
println(@test DP.normsq(dX) == DP.normsq(dX.rr) )

##

@info("performance/allocation test ")

function bm_copy!(Y, X, a)
   for i = 1:length(Y)
      Y[i] = X[i] * a 
   end 
   return Y
end

PositionStateF64 = typeof(PState(rr = zero(SVector{3, Float64})))

Xs = [ rand(PositionStateF64) for _=1:100 ]
Ys = [ zero(PositionStateF64) for _=1:100 ]
a = rand() 
let Ys = Ys, Xs = Xs, a = a
   bm_copy!(Ys, Xs, a)
   nalloc = @allocated bm_copy!(Ys, Xs, a)
   @show nalloc 
   println(@test nalloc == 0) 
end

