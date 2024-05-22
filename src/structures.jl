
import AtomsBase 
import AtomsBase: AbstractSystem, ChemicalElement, 
                  position, velocity, atomic_mass, atomic_number, 
                  atomic_symbol

# --------------------------------------------------- 
# an `Atom` is now just a `PState`, so we define 
# accessors for the PState fields with canonical names. 

symbol(::typeof(position)) = :𝐫 
symbol(::typeof(velocity)) = :𝐯
symbol(::typeof(atomic_mass)) = :m
symbol(::typeof(atomic_symbol)) = :Z

position(atom::PState) = atom.𝐫 
velocity(atom::PState) = atom.𝐯
atomic_mass(atom::PState) = atom.m
atomic_symbol(atom::PState) = atom.Z   # this one I'm not sure about 

atomic_number(atom::PState) = atomic_number(atomic_symbol(atom))

""" 
Generate an atom with the given properties. 
"""
atom(at; properties = (position, atomic_mass, atomic_symbol)) = 
      PState((; [symbol(p) => p(at) for p in properties]...))

# --------------------------------------------------- 

mutable struct AosSystem{D, TCELL, TPART} <: AbstractSystem{D} 
   cell::TCELL 
   particles::Vector{TPART}
   # -------- 
   meta::Dict{String, Any}
end


function AosSystem(sys::AbstractSystem; 
                   properties = (position, atomic_mass, atomic_symbol), )

   X = [ atom(sys[i]; properties = properties) for i = 1:length(sys) ]
   cell = AtomsBase.get_cell(sys)
   D = AtomsBase.n_dimensions(cell)
   return AosSystem{D, typeof(cell), eltype(X)}(cell, X, Dict{String, Any}())
end


# ---------------------------------------------------
# implementing the interface 

Base.length(at::AosSystem) = length(at.particles)

Base.getindex(at::AosSystem, i::Int) = at.particles[i]

for f in (:position, :velocity, :atomic_mass, :atomic_symbol)
   @eval $f(sys::AosSystem) = [ $f(x) for x in sys.particles ]
   @eval $f(sys::AosSystem, i::Integer) = $f(sys.particles[i]) 
   @eval $f(sys::AosSystem, inds::AbstractVector) = [$f(sys.particles[i]) for i in inds]
end

AtomsBase.get_cell(at::AosSystem) = at.cell

AtomsBase.n_dimensions(at::AosSystem) = AtomsBase.n_dimensions(at.cell)
AtomsBase.bounding_box(at::AosSystem) = AtomsBase.bounding_box(at.cell)
AtomsBase.boundary_conditions(at::AosSystem) = AtomsBase.boundary_conditions(at.cell)
AtomsBase.periodicity(at::AosSystem) = AtomsBase.periodicity(at.cell)

