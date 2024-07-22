
import AtomsBase 
import AtomsBase: AbstractSystem,  
                  position, velocity, atomic_mass, atomic_number, 
                  atomic_symbol, n_dimensions, bounding_box, 
                  boundary_conditions, periodicity

import DecoratedParticles.Tmp: SystemWithCell, PCell, ChemicalElement, get_cell 

export get_cell, AosSystem, SoaSystem, ChemicalElement


# --------------------------------------------------- 
# an `Atom` is now just a `PState`, so we define 
# accessors for the PState fields with canonical names. 

# list of properties that DP knows about: 
_list_of_properties = [ 
   (:𝐫, :position),
   (:𝐯, :velocity),
   (:𝑚, :atomic_mass),
   (:𝑍, :atomic_symbol),
   (:𝑞, :charge),
   (:μ, :dipole), 
   (:𝐩, :momentum),
   (:𝐸, :energy),
   (:𝑀, :mass),
   ]

for (sym, name) in _list_of_properties
   @eval $name(x::XState) = x.$sym
   @eval symbol(::typeof($name)) = $(Meta.quot(sym))
end

_post(p, a) = a 
_post(::typeof(atomic_symbol), sym::Symbol) = ChemicalElement(sym)
_post(::typeof(atomic_symbol), z::Integer) = ChemicalElement(z)

""" 
Generate an atom with the given properties. 
"""
atom(at; properties = (position, atomic_mass, atomic_symbol)) = 
      PState((; [symbol(p) => _post(p, p(at)) for p in properties]...))

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

# for (sym, name) in [ _list_of_properties; [(:ignore, :atomic_number) ] ]
for (sym, name) in _list_of_properties
   @eval $name(sys::AosSystem) = [ $name(x) for x in sys.particles ]
   @eval $name(sys::AosSystem, i::Integer) = $name(sys.particles[i])
   @eval $name(sys::AosSystem, inds::AbstractVector) = [$name(sys.particles[i]) for i in inds]
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

   arrays = (; [symbol(p) => _post.(p, p(sys)) for p in properties]... )
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

# for (sym, name) in [ _list_of_properties; [(:ignore, :atomic_number) ] ]
for (sym, name) in _list_of_properties
   @eval $name(sys::SoaSystem) = copy(sys.arrays.$sym)
   @eval $name(sys::SoaSystem, i::Integer) = sys.arrays.$sym[i]
   @eval $name(sys::SoaSystem, inds::AbstractVector) = sys.arrays.$sym[inds]
end

get_cell(at::SoaSystem) = at.cell

for f in (:n_dimensions, :bounding_box, :boundary_conditions, :periodicity)
   @eval $f(at::SoaSystem) = $f(at.cell)
end


# ---------------------------------------------------
# atomic_number overload 

atomic_number(atom::PState) = 
      atomic_number(atomic_symbol(atom))

atomic_number(sys::Union{SoaSystem, AosSystem}, i::Integer) = 
      atomic_number(atomic_symbol(sys, i))

atomic_number(sys::Union{SoaSystem, AosSystem}, args...) = 
      atomic_number.(atomic_symbol(sys, args...))



# ---------------------------------------------------------------
#  Extension of the AtomsBase interface with setter functions 

for (sym, name) in _list_of_properties
   set_name = Symbol("set_$name")
   set_name_ip = Symbol("set_$(name)!")
   set_names_ip = Symbol("set_$(name)s!")
   @eval export $set_name, $set_name_ip, $set_names_ip 
   # e.g. set_position(x::PState, 𝐫) = setproperty(x, :𝐫, 𝐫)
   @eval $set_name(x::PState, t) = setproperty(x, $(Meta.quot(sym)), t)
   # e.g. set_position!(sys, i, 𝐫)
   @eval begin
      function $set_name_ip(sys::AosSystem, i::Integer, t) 
         xi = sys.particles[i] 
         sys.particles[i] = $set_name(xi, t)
         return nothing 
      end
   end
   @eval begin 
      function $set_name_ip(sys::SoaSystem, i::Integer, t)
         sys.arrays.$sym[i] = t
         return nothing 
      end
   end
   # e.g. set_positions!(sys, X)
   @eval begin
      function $set_names_ip(sys::AosSystem, Rs) 
         for i = 1:length(sys) 
            $set_name_ip(sys, i, Rs[i])
         end
         return sys  
      end
   end
   @eval begin 
      function $set_names_ip(sys::SoaSystem, R::AbstractVector)
         copy!(sys.arrays.$sym, R)
         return sys  
      end
   end
end

# todo: set_atomic_number & co 

# cell-level setter 

export set_bounding_box!

function set_bounding_box!(sys::Union{SoaSystem, AosSystem}, bb)
   cell = PCell(bb, periodicity(sys))
   sys.cell = cell
   return nothing 
end