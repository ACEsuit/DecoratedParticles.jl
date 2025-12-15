# For interactive debugging only 
# using Pkg; Pkg.activate(joinpath(@__DIR__(), ".."))
# using TestEnv; TestEnv.activate();
# Pkg.develop(url = joinpath(@__DIR__(), ".."))

##

@info("AtomsBase tests") 
using DecoratedParticles, AtomsBase, StaticArrays, Unitful, Test, AtomsBuilder 
using AtomsBase: Atom, ChemicalSpecies, cell_vectors, periodicity
import DecoratedParticles as DP 

using LinearAlgebra: I 

## 
#generate an atom and check that the accessors work

Z0 = 6 
x = PState(𝐫 = SA[1.0, 2.0, 3.0], 𝐯 = SA[0.1, 0.2, 0.3], 
           𝑚 = 1.0, S = ChemicalSpecies(Z0) )
display(x)           
@test position(x) == x.𝐫
@test velocity(x) == x.𝐯
@test mass(x) == x.𝑚
# @test atomic_symbol(x) == x.S
@test atomic_number(x) == Z0


## 
#convert an Atom 

at = Atom(6, SA[1.0, 2.0, 3.0]u"Å"; mass = 1.0u"u")
x = DP.atom(at; properties = (position, mass, species))
display(x)
@test x.𝐫 == position(x) == position(at)
@test x.𝑚 == mass(x) == mass(at)
@test x.S == species(x)
@test DP.symbol(position) == :𝐫
@test DP.symbol(mass) == :𝑚
@test DP.symbol(atomic_number) == :𝑍
@test DP.symbol(species) == :S

## 
# convert an entire system 

sys = rattle!(bulk(:Si, cubic=true) * 2, 0.1);
aos = DP.AosSystem(sys)
soa = DP.SoaSystem(sys)

aos[1]
soa[1]

for i = 1:10 
   @test aos[i] == soa[i]

   for f in (position, mass, species)
      @test f(aos, i) == f(soa, i)
   end 
end 

for f in (cell, periodicity, cell_vectors, n_dimensions)
   @test f(aos) == f(soa)
end

s14 = ChemicalSpecies(14)
for _sys in (aos, soa)
   @test species(_sys, :) == fill(s14, length(_sys))
   @test species(_sys, 5) == s14 
   @test species(_sys, [2,4,7]) == fill(s14, 3)
end

## 
# some performance related tests 

sys = rattle!(bulk(:Si, cubic=true) * 2, 0.1);
aos = DP.AosSystem(sys)
soa = DP.SoaSystem(sys)

x1 = aos[1] 
x2 = soa[1]
@test isbits(x1)           
@test isbits(x2)

@info("Checking allocations during accessors")
_check_allocs(sys) = ( (@allocated position(sys, 1)) + 
                       (@allocated mass(sys, 1) ) +
                       (@allocated sys[1] ) )
@test _check_allocs(aos) == 0
@test _check_allocs(soa) == 0 


## 
# setters 

sys = rattle!(bulk(:Si, cubic=true) * 2, 0.1);
aos = DP.AosSystem(sys)

x = aos[1]
𝐫 = x.𝐫
𝐫1 = 1.01 * 𝐫
x1 = DP.setproperty(x, :𝐫, 𝐫1)
@test x1.𝐫 == 𝐫1
@test x1.𝑚 == x.𝑚
@test x1.S == x.S

x2 = DP.set_position(x, 𝐫1)
@test x2 == x1

##

X = position(sys, :) + 0.01 * randn(SVector{3, Float64}, length(aos)) * u"Å"
aos1 = deepcopy(aos)
DP.set_positions!(aos1, X)
@test all(position(aos1, :) .== X)
@test !any(position(aos, :) .== X)

DP.set_position!(aos1, 1, 𝐫1)
@test position(aos1, 1) == 𝐫1

@allocated DecoratedParticles.set_position(x, 𝐫1)
@test (@allocated DP.set_positions!(aos1, X)) == 0

## 

soa = DP.SoaSystem(aos)
@test all(position(soa, :) .== position(aos, :))
DP.set_positions!(soa, X)
@test all(position(soa, :) .== position(aos1, :))
DP.set_position!(soa, 1, 𝐫1)
@test position(soa, 1) == 𝐫1

@test (@allocated DP.set_positions!(soa, X)) == 0 

## 

bb = cell_vectors(soa)
bb1 = ntuple(i -> (I + 0.01*randn(SMatrix{3,3,Float64})) * bb[i], 3)
DP.set_cell_vectors!(soa, bb1)
@test cell_vectors(soa) == bb1
DP.set_cell_vectors!(aos, bb1)
@test cell_vectors(aos) == bb1

