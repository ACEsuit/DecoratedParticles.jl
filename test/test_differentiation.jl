
using DecoratedParticles, StaticArrays, Test, ForwardDiff
using DecoratedParticles: PState, VState, grad_fd

import DecoratedParticles as DP

##

@info("Testing differentiation tooling")
@info(" ... preliminaries")

rand_nt() = (q = randn(), r = randn(SVector{3, Float64}), z = rand(1:10))

x = rand_nt()
v_nt = DP._ctsnt(x)
@test v_nt == (q = x.q, r = x.r)
v = DP._nt2svec(v_nt)
@test v == SVector{4, Float64}(x.q, x.r...)
v_nt1 = DP._svec2nt(v, v_nt)
@test v_nt1 == v_nt

##

@info(" ... grad_fd")

module TestDiff
   using StaticArrays, ForwardDiff
   import DecoratedParticles: PState, VState

   struct F{N, T}; W::SVector{N, T}; end

   # random expression, but representative in terms of simplicity
   (f::F)(x) = sum(x.r .* x.r) * x.q / (1 + f.W[x.z]^2)

   # manual gradient, as NamedTuple
   grad_man(f::F, x) = ( r2 = sum(x.r .* x.r); w = 1 / (1 + f.W[x.z]^2);
                        (q = r2 * w, r = 2 * x.r * x.q * w, )  )

   # gradient via ForwardDiff, component by component
   function grad_1(f::F, x)
      ∂q = ForwardDiff.derivative(q -> f((q=q,   r=x.r, z=x.z)), x.q)
      ∂r = ForwardDiff.gradient(r -> f((q=x.q, r=r,   z=x.z)), x.r)
      return (q = ∂q, r = ∂r)
   end
end

f = TestDiff.F(@SVector randn(10))

for ntest = 1:20
   local x, x_dp
   x = rand_nt()
   x_dp = PState(; x...)
   g_man = TestDiff.grad_man(f, x)
   g_fd1 = TestDiff.grad_1(f, x)
   g_nt = grad_fd(f, x)          # NamedTuple input -> NamedTuple output
   g_dp = grad_fd(f, x_dp)       # PState input -> VState output
   @test all(g_fd1[sym] ≈ g_man[sym] for sym in fieldnames(typeof(g_man)))
   @test all(g_nt[sym] ≈ g_man[sym] for sym in fieldnames(typeof(g_man)))
   @test g_dp isa VState
   @test g_dp ≈ VState(g_man)
end

##

@info(" ... NamedTuple path vs PState path consistency")

# the NamedTuple path (_ctsnt/_iscts) and the XState path
# (vstate_type/_findcts) must agree on which fields are continuous,
# including Complex and Dual fields
using ForwardDiff: Dual

nt_mixed = (q = randn(), r = randn(SVector{3, Float64}), z = rand(1:10),
            c = randn(ComplexF64), d = Dual(randn(), 1.0))
@test keys(DP._ctsnt(nt_mixed)) == DP._ctssyms(PState(; nt_mixed...))

# Dual-valued fields survive the PState -> VState conversion
# (needed for nested ForwardDiff through states)
x_dual = PState(q = Dual(randn(), 1.0), r = randn(SVector{3, Float64}),
                z = rand(1:10))
@test haskey(getfield(VState(x_dual), :x), :q)

##

@info(" ... allocation tests")

N = 1000
X = [ rand_nt() for _ in 1:N ]
TG = typeof(grad_fd(f, X[1]))
gY = Vector{TG}(undef, N)

function count_alloc(f, gY, X)
   gfun = x -> grad_fd(f, x)
   map!(gfun, gY, X)
   return @allocations map!(gfun, gY, X)
end

@test count_alloc(f, gY, X) <= 1
