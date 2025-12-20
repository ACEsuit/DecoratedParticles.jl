

using DecoratedParticles, StaticArrays, Test, LinearAlgebra
using DecoratedParticles: PState, VState
import DecoratedParticles as DP 

##

X = PState( rr = randn(SVector{3, Float64}) )

SYMS = DP._syms(X)
println(@test SYMS == (:rr,))

TT = DP._tt(X)
println(@test( TT == Tuple{SVector{3, Float64}} ))

TDX = DP.vstate_type(X)
println(@test TDX == VState{NamedTuple{(:rr,), Tuple{SVector{3, Float64}}}})
println(@test TDX == DP.vstate_type(typeof(X)))

dX = VState(X)
println(@test typeof(dX) == DP.vstate_type(X))
println(@test typeof(dX) == DP.vstate_type(typeof(dX)))
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

## 

x = PState(a = 1.0, b = 2, c = true)
y = DP.setproperty(x, :a, 2.3)
z = DP.setproperty(x, Val(:a), 2.3)
@test y.a == 2.3 && y.b == x.b && y.c == x.c
@test z.a == 2.3 && z.b == x.b && z.c == x.c
@test y == z

function transform_states!(Xs, sym)
   for i = 1:length(Xs)
      xi = Xs[i] 
      Xs[i] = DecoratedParticles.setproperty(xi, sym, rand())
   end
   return nothing 
end

# test that there are no allocations 
Xs = [ PState(a = 1.0, b = 2, c = true)  for _ = 1:1000 ]
transform_states!(Xs, :a)
@test (@allocated transform_states!(Xs, :a)) == 0
transform_states!(Xs, Val{:a}())
@test (@allocated transform_states!(Xs, Val{:a}())) == 0 

# the two codes are essentially equivalent!!
# using BenchmarkTools
# @btime transform_states!($Xs, $(:a))
# @btime transform_states!($Xs, Val{:a}())

