


# --------------------------------------------- 
#  Simple wrapper for chemical element type 

import PeriodicTable
using Unitful

import Base: ==, convert, show, length
import AtomsBase: atomic_mass, atomic_symbol, atomic_number, element 

"""
Encodes a chemical element by wrapping an integer that represents the atomic 
number. 
"""
struct ChemicalElement
    atomic_number::UInt8
end

Base.show(io::IO, element::ChemicalElement) = 
      print(io, Symbol(element))

ChemicalElement(sym::Symbol) = ChemicalElement(_sym2z[sym])  
ChemicalElement(sym::ChemicalElement) = sym

==(a::ChemicalElement, sym::Symbol) = (Symbol(a) == sym)

# -------- fast access to the periodic table 

@assert length(PeriodicTable.elements) == maximum(el.number for el in PeriodicTable.elements)
@assert all(el.number == i for (i, el) in enumerate(PeriodicTable.elements))

const _chem_el_info = [ 
      (symbol = Symbol(PeriodicTable.elements[z].symbol), 
       atomic_mass = PeriodicTable.elements[z].atomic_mass, ) 
       for z in 1:length(PeriodicTable.elements)
      ]

const _sym2z = Dict{Symbol, UInt8}() 
for z in 1:length(_chem_el_info)
   _sym2z[_chem_el_info[z].symbol] = z
end

# -------- accessor functions 

# UInt8 is not readable
atomic_number(element::ChemicalElement) = Int(element.atomic_number)  

atomic_symbol(element::ChemicalElement) = element 

Base.convert(::Type{Symbol}, element::ChemicalElement) = Symbol(element) 
Base.Symbol(element::ChemicalElement) = _chem_el_info[element.atomic_number].symbol

atomic_mass(element::ChemicalElement) = _chem_el_info[element.atomic_number].atomic_mass

rich_info(element::ChemicalElement) = PeriodicTable.elements[element.atomic_number]

element(element::ChemicalElement) = rich_info(element)

