# DecoratedParticles.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ACEsuit.github.io/DecoratedParticles.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ACEsuit.github.io/DecoratedParticles.jl/dev/)
[![Build Status](https://github.com/ACEsuit/DecoratedParticles.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ACEsuit/DecoratedParticles.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is a small package, spun out of [`ACE.jl`](https://github.com/ACEsuit/ACE.jl). The original intended use case is managing (lists of) decorated particles, i.e., point clouds embedded in some vector space, where each point is decorated with additional features such as chemical species, charge, mass, etc. 

Example usage:
```julia
using DecoratedParticles, StaticArrays
# a silicon atom
x = PState(𝐫 = randn(SVector{3, Float64}), z = 14)
# ⟨𝐫:[-0.91, -0.87, -0.42], z:14⟩
# extract its position
x.𝐫
# 3-element SVector{3, Float64} with indices SOneTo(3):
#  -1.1469536186585183
#  -0.1832512302259138
#   1.0216715637205427
```
