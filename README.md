# DecoratedParticles.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ACEsuit.github.io/DecoratedParticles.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ACEsuit.github.io/DecoratedParticles.jl/dev/) -->
[![Build Status](https://github.com/ACEsuit/DecoratedParticles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ACEsuit/DecoratedParticles.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a small package, spun out of [`ACE.jl`](https://github.com/ACEsuit/ACE.jl). The original intended use case is managing (lists of) decorated particles, i.e., point clouds embedded in some vector space, where each point is decorated with additional features such as chemical species, charge, mass, etc. Documentation for now is this readme. 

### Example usage

#### `PState` and `VState` 

```julia
using DecoratedParticles, StaticArrays, LinearAlgebra, Zygote 
using DecoratedParticles: PState, VState
import DecoratedParticles as DP 

x1 = PState( 𝐫 = randn(SVector{3, Float64}), z = 14 )
# 〖𝐫:[-0.74, -2.27, -0.83], z:14〗
x2 = PState( 𝐫 = randn(SVector{3, Float64}), z = 14 )
# 〖𝐫:[-0.63, 0.67, -0.56], z:14〗
𝐫12 = VState(x2 - x1)
# ｟𝐫:[0.11, 2.94, 0.27]｠

# extract the position 
x1.𝐫
# 3-element SVector{3, Float64} with indices SOneTo(3):
#  -0.7424735839283951
#  -2.271376247109223
#  -0.8265064008465374

# arithmetic on particle states 
x1 + 𝐫12 ≈ x2
# true 

f(X) = sum(DP.normsq(x.𝐫) for x in X)
f([x1, x2])
# 4.115...

# the gradient of a PState is a VState 
g = Zygote.gradient(f, [x1, x2])[1] 
# 2-element Vector{VState{@NamedTuple{𝐫::SVector{3, Float64}}}}:
#｟𝐫:[-1.48, -4.54, -1.65]｠
#｟𝐫:[-1.26, 1.35, -1.12]｠

g[1].𝐫 ≈ 2 * x1.𝐫
# true 

# Some property symbols are standardized, e.g. 𝐫 always means `position`
x1.𝐫 == position(x1)   # true 

# a 4-momentum might look like this 
p = PState(𝐩 = randn(SVector{3, Float64}), 𝐸 = rand()) 
p.𝐩 == DP.momentum(p)
p.𝐸 == DP.energy(p)
```

### Prototype AtomsBase system implementations 

AtomsBase system implementations are provided via a package extension. 
Both AosSystem and SoaSystem are fully flexible regarding the 
properties of the particles but have the same performance as `FastSystem`.

```julia
using AtomsBase, AtomsBuilder
sys = rattle!(bulk(:Si, cubic=true) * 2, 0.1);   # AtomsBase.FlexibleSystem
fsys = FastSystem(sys);                          # AtomsBase.FastSystem
aos = DP.aos_system(sys);
soa = DP.soa_system(sys);

x1 = aos[1]   # PState, just sys.particles[1]
x2 = soa[1]   # PState, generated from the arrays in sys 
isbits(x1)    # true 
isbits(x2)    # true 

display(x1)   # 〖𝐫:[-0.01, -0.02, -0.1] Å, 𝑚:28.085 u, 𝑍:Si〗

# specific symbols are taken equivalent to AtomsBase accessor functions e.g. 
position(x1) == x1.𝐫       # true 
mass(x1) == x1.𝑚    # true
species(x1) == x1.S  # true

# Performance

# accessors are non-allocating: 
_check_allocs(sys) =  ( (@allocated position(sys, 1)) + 
                        (@allocated mass(sys, 1)) + 
                        (@allocated sys[1])  )
_check_allocs(sys)   # 288 
_check_allocs(fsys)  # 0
_check_allocs(aos)   # 0 
_check_allocs(soa)   # 0 

# this has performance implications
using BenchmarkTools

# Silly test 1 : sum up the positions via `position(sys, i)` accessor
silly_test_1(sys) = sum( position(sys, i) for i = 1:length(sys) )
@btime silly_test_1($sys)   #   8.819 μs (320 allocations: 12.00 KiB)
@btime silly_test_1($fsys)  #   50.405 ns (0 allocations: 0 bytes)
@btime silly_test_1($aos)   #   50.447 ns (0 allocations: 0 bytes)
@btime silly_test_1($soa)   #   50.405 ns (0 allocations: 0 bytes)

silly_test_2(sys) = sum( position(x) for x in sys )
@btime silly_test_2($sys)   #   10.750 μs (256 allocations: 18.00 KiB)
@btime silly_test_2($fsys)  #   48.118 ns (0 allocations: 0 bytes)
@btime silly_test_2($aos)   #   47.950 ns (0 allocations: 0 bytes)
@btime silly_test_2($soa)   #   48.794 ns (0 allocations: 0 bytes)
```
