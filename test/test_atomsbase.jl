
using DecoratedParticles, AtomsBase, StaticArrays, Unitful, Test 
using AtomsBase: ChemicalElement, Atom 
using AtomsBuilder: bulk, rattle! 
DP = DecoratedParticles

## 
#generate an atom and check that the accessors work

x = PState(𝐫 = SA[1.0, 2.0, 3.0], 𝐯 = SA[0.1, 0.2, 0.3], 
           m = 1.0, Z = ChemicalElement(6) )
display(x)           
@test position(x) == x.𝐫
@test velocity(x) == x.𝐯
@test atomic_mass(x) == x.m
@test atomic_symbol(x) == x.Z
@test atomic_number(x) == 6

## 
#convert an Atom 

at = Atom(6, SA[1.0, 2.0, 3.0]u"Å"; atomic_mass = 1.0u"u")
x = DP.atom(at; properties = (position, atomic_mass, atomic_symbol))
display(x)
@test x.𝐫 == position(x) == position(at)
@test x.m == atomic_mass(x) == atomic_mass(at)
@test x.Z == atomic_symbol(x) == atomic_symbol(at)


## 
# convert an entire system 

sys = rattle!(bulk(:Si, cubic=true) * 2, 0.1);
aos = DP.AosSystem(sys);
soa = DP.SoaSystem(sys);

for i = 1:10 
   @test aos[i] == soa[i]

   for f in (position, atomic_mass, atomic_symbol)
      @test f(aos, i) == f(soa, i)
   end 
end 

for f in (get_cell, periodicity, boundary_conditions, bounding_box, n_dimensions)
   @test f(aos) == f(soa)
end

