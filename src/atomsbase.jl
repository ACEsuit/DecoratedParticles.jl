
import AtomsBase 
import AtomsBase: AbstractSystem,  
                  position, velocity, atomic_mass, atomic_number, 
                  atomic_symbol, n_dimensions, bounding_box, 
                  boundary_conditions, periodicity

import DecoratedParticles.Tmp: SystemWithCell, PCell, ChemicalElement, get_cell 

# --------------------------------------------------- 
# an `Atom` is now just a `PState`, so we define 
# accessors for the PState fields with canonical names. 

symbol(::typeof(position)) = :𝐫 
symbol(::typeof(velocity)) = :𝐯
symbol(::typeof(atomic_mass)) = :𝑚
symbol(::typeof(atomic_symbol)) = :𝑍

const _atom_syms = (𝐫 = position, 
                    𝐯 = velocity, 
                    𝑚 = atomic_mass, 
                    𝑍 = atomic_symbol)

position(atom::PState) = atom.𝐫 
velocity(atom::PState) = atom.𝐯
atomic_mass(atom::PState) = atom.𝑚
atomic_symbol(atom::PState) = atom.𝑍   # this one I'm not sure about 

atomic_number(atom::PState) = atomic_number(atomic_symbol(atom))

_post(a) = a 
_post(sym::Symbol) = ChemicalElement(sym)

""" 
Generate an atom with the given properties. 
"""
atom(at; properties = (position, atomic_mass, atomic_symbol)) = 
      PState((; [symbol(p) => _post(p(at)) for p in properties]...))

# --------------------------------------------------- 
# Array of Structs System 
#
mutable struct AosSystem{D, TCELL, TPART} <: SystemWithCell{D, TCELL} 
   cell::TCELL 
   particles::Vector{TPART}
   # -------- 
   meta::Dict{String, Any}
end


function AosSystem(sys::AbstractSystem; 
                   properties = (position, atomic_mass, atomic_symbol), )

   X = [ atom(sys[i]; properties = properties) for i = 1:length(sys) ]
   cell = get_cell(sys)
   D = AtomsBase.n_dimensions(cell)
   return AosSystem{D, typeof(cell), eltype(X)}(cell, X, Dict{String, Any}())
end


# implementing the AtomsBase interface 

Base.length(at::AosSystem) = length(at.particles)

Base.getindex(at::AosSystem, i::Int) = at.particles[i]
Base.getindex(at::AosSystem, inds::AbstractVector) = at.particles[inds]

for f in (:position, :velocity, :atomic_mass, :atomic_symbol)
   @eval $f(sys::AosSystem) = [ $f(x) for x in sys.particles ]
   @eval $f(sys::AosSystem, i::Integer) = $f(sys.particles[i]) 
   @eval $f(sys::AosSystem, inds::AbstractVector) = [$f(sys.particles[i]) for i in inds]
end

# AtomsBase. 
get_cell(at::AosSystem) = at.cell

for f in (:n_dimensions, :bounding_box, :boundary_conditions, :periodicity)
   @eval $f(at::AosSystem) = $f(at.cell)
end


# --------------------------------------------------- 
# Struct of Arrays System 

mutable struct SoaSystem{D, TCELL, NT} <: SystemWithCell{D, TCELL} 
   cell::TCELL 
   arrays::NT 
   # -------- 
   meta::Dict{String, Any}
end

function SoaSystem(sys::AbstractSystem; 
                   properties = (position, atomic_mass, atomic_symbol), )

   arrays = (; [symbol(p) => _post.(p(sys)) for p in properties]... )
   cell = get_cell(sys)
   D = AtomsBase.n_dimensions(cell)
   return SoaSystem{D, typeof(cell), typeof(arrays)}(
                  cell, arrays, Dict{String, Any}())
end


# implementing the AtomsBase interface

Base.length(at::SoaSystem) = length(at.arrays[1])

# this implementation seems canonical but appears to be type unstable 
# unclear to me exactly why, maybe it can be fixed. 
# function Base.getindex(sys::SoaSystem{D, TCELL, NT}, i::Integer) where {D, TCELL, NT}
#    SYMS = _syms(NT)
#    return PState(; ntuple(a -> SYMS[a] => sys.arrays[SYMS[a]][i], length(SYMS))...)
# end

@generated function Base.getindex(sys::SoaSystem{D, TCELL, NT}, i::Integer) where {D, TCELL, NT}
   SYMS = _syms(NT) 
   # very naive code-writing ... probably there is a nicer way ... 
   code = "PState("
   for sym in SYMS 
      code *= "$(sym) = sys.arrays.$sym[i], " 
   end
   code *= ")"
   return quote 
      $(Meta.parse(code))
   end
end

Base.getindex(sys::SoaSystem, inds::AbstractVector{<: Integer}) = 
      [ sys[i] for i in inds ]

for f in (:position, :velocity, :atomic_mass, :atomic_symbol)
   @eval $f(sys::SoaSystem) = getfield(sys.arrays, symbol($f))
   @eval $f(sys::SoaSystem, i::Integer) = getfield(sys.arrays, symbol($f))[i]
   @eval $f(sys::SoaSystem, inds::AbstractVector) = getfield(sys.arrays, symbol($f))[inds]
end

# AtomsBase.
get_cell(at::SoaSystem) = at.cell

for f in (:n_dimensions, :bounding_box, :boundary_conditions, :periodicity)
   @eval $f(at::SoaSystem) = $f(at.cell)
end




# ---------------------------------------------------------------
#  Extension of the AtomsBase interface with setter functions 

