using DecoratedParticles, StaticArrays, LinearAlgebra, Zygote 
using DecoratedParticles: PState, VState
DP = DecoratedParticles

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
