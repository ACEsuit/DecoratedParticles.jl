
using DecoratedParticles, AtomsBase, StaticArrays, Unitful, Test, AtomsBuilder 
using AtomsBase: Atom 
using DecoratedParticles.Tmp: ChemicalElement, get_cell 
DP = DecoratedParticles

using LinearAlgebra: I 

## 
#generate an atom and check that the accessors work

x = PState(𝐫 = SA[1.0, 2.0, 3.0], 𝐯 = SA[0.1, 0.2, 0.3], 
           𝑚 = 1.0, 𝑍 = ChemicalElement(6) )
display(x)           
@test position(x) == x.𝐫
@test velocity(x) == x.𝐯
@test atomic_mass(x) == x.𝑚
@test atomic_symbol(x) == x.𝑍
@test atomic_number(x) == 6

## 
#convert an Atom 

at = Atom(6, SA[1.0, 2.0, 3.0]u"Å"; atomic_mass = 1.0u"u")
x = DP.atom(at; properties = (position, atomic_mass, atomic_symbol))
display(x)
@test x.𝐫 == position(x) == position(at)
@test x.𝑚 == atomic_mass(x) == atomic_mass(at)
@test x.𝑍 == atomic_symbol(x) == atomic_symbol(at)


## 
# convert an entire system 

sys = rattle!(bulk(:Si, cubic=true) * 2, 0.1);
aos = DP.AosSystem(sys)
soa = DP.SoaSystem(sys)

aos[1]
soa[1]

for i = 1:10 
   @test aos[i] == soa[i]

   for f in (position, atomic_mass, atomic_symbol)
      @test f(aos, i) == f(soa, i)
   end 
end 

for f in (get_cell, periodicity, boundary_conditions, bounding_box, n_dimensions)
   @test f(aos) == f(soa)
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
                       (@allocated atomic_mass(sys, 1) ) +
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
x1 = DP.set_property(x, :𝐫, 𝐫1)
@test x1.𝐫 == 𝐫1
@test x1.𝑚 == x.𝑚
@test x1.𝑍 == x.𝑍

x2 = DP.set_position(x, 𝐫1)
@test x2 == x1

##

X = position(sys) + 0.01 * randn(SVector{3, Float64}, length(aos)) * u"Å"
aos1 = deepcopy(aos)
DP.set_positions!(aos1, X)
@test all(position(aos1) .== X)
@test !any(position(aos) .== X)

DP.set_position!(aos1, 1, 𝐫1)
@test position(aos1, 1) == 𝐫1

@allocated DecoratedParticles.set_position(x, 𝐫1)
@test (@allocated DP.set_positions!(aos1, X)) == 0

## 

soa = DP.SoaSystem(aos)
@test all(position(soa) .== position(aos))
DP.set_positions!(soa, X)
@test all(position(soa) .== position(aos1))
DP.set_position!(soa, 1, 𝐫1)
@test position(soa, 1) == 𝐫1

@test (@allocated DP.set_positions!(soa, X)) == 0 

## 

bb = bounding_box(soa)
bb1 = ntuple(i -> (I + 0.01*randn(SMatrix{3,3,Float64})) * bb[i], 3)
DP.set_bounding_box!(soa, bb1)
@test bounding_box(soa) == bb1
DP.set_bounding_box!(aos, bb1)
@test bounding_box(aos) == bb1

