

export PState, VState

using StaticArrays, NamedTupleTools

import Base: *, +, -, zero, rand, randn, show, promote_rule, rtoldefault, 
             isapprox, getproperty, real 

import LinearAlgebra: norm, promote_leaf_eltypes

export PState, VState, vstate_type, setproperty 

abstract type XState{NT <: NamedTuple} end 

"""
`struct State` the main type for states of input variables (particles). 
This type is intended only for storing of information but no arithmetic 
should be performed on it. For the latter, we have the DState.
"""
struct PState{NT <: NamedTuple} <: XState{NT}
   x::NT

   # standard constructor - SYMS are the same 
   PState{NT}(t::NT1) where {NT <: NamedTuple{SYMS}, NT1 <: NamedTuple{SYMS}} where {SYMS} = 
         new{NT1}(t)

   # if SYMS are not the same we automatically merge them
   PState{NT}(t::NT1) where {NT <: NamedTuple, NT1 <: NamedTuple} = 
         PState( merge( _x(zero(State{NT})), t ) )
end

"""
`struct VState`: A `State`-like variable but acting like a vector with arithmetic 
operations defined on it, while `State` acts more like a Point i.e. a fixed object that cannot 
be manipulated. The main application of `VState` is as a derivative of a 
`State`; see also `Vstate_type`, ... 
"""
struct VState{NT <: NamedTuple} <: XState{NT}
   x::NT

   # standard constructor - SYMS are the same 
   VState{NT}(t::NT1) where {NT  <: NamedTuple{SYMS}, 
                             NT1 <: NamedTuple{SYMS}} where {SYMS} = 
         new{NT1}(t)

   # if SYMS are not the same we automatically merge them
   VState{NT}(t::NT1) where {NT <: NamedTuple, NT1 <: NamedTuple} = 
         VState( merge( _x(zero(VState{NT})), t ) )

   VState{NT}(dX::VState) where {NT} = VState{NT}(_x(dX))
end

# ------------------------------------------ 
#  Further PState Constructors 

# the two standard outward facing constructors
PState(nt::NT) where {NT <: NamedTuple} = PState{NT}(nt)

PState(; kwargs...) = PState(NamedTuple(kwargs))

# accessing X.nt and the fields of X.nt 
# this relies heavily on constant propagation 
_x(X::XState) = getfield(X, :x)

getproperty(X::XState, sym::Symbol) = getproperty(_x(X), sym)

Base.length(::XState) = 1 


# ----------- some basic manipulations 

_syms(nt::NamedTuple{SYMS, TT}) where {SYMS, TT} = SYMS
_syms(nt::Type{NamedTuple{SYMS, TT}}) where {SYMS, TT} = SYMS


# extract the symbols and the types 
_syms(X::XState) = _syms(typeof(X))
_syms(::Type{<: XState{NamedTuple{SYMS, TT}}}) where {SYMS, TT} = SYMS

@generated function _syms(x1::TX1, x2::TX2) where {TX1 <: XState, TX2 <: XState}
   vsyms = unique( tuple(_syms(TX1)..., _syms(TX2)...) )
   syms = tuple(vsyms...)
   quote
      $syms
   end
end



_tt(X::XState) = _tt(typeof(X))
_tt(::Type{<: XState{NamedTuple{SYMS, TT}}}) where {SYMS, TT} = TT

_symstt(X::XState) = _symstt(typeof(X))
_symstt(::Type{<: XState{NamedTuple{SYMS, TT}}}) where {SYMS, TT} = SYMS, TT 


# --------------------------------------------------
#   specify cts vs categorical variables 

# which properties are continuous - this is important since certain operations 
# act only on continuous variables but not on categorical variables. 

const CTSTT = Union{AbstractFloat, 
                    Complex{<: AbstractFloat},
                    SVector{N, <: AbstractFloat}, 
                    SVector{N, <: Complex}} where {N}

"""
Find the indices of continuous properties, and return the 
indices as well as the symbols and types
"""
_findcts(X::TX) where {TX <: XState} = _findcts(TX)

@generated function _findcts(::Type{TX}) where {TX <: XState}
   SYMS, TT = _symstt(TX)
   icts = findall(T -> T <: CTSTT, TT.types)
   quote 
      $(SVector(icts...))
   end
end

_ctssyms(X::TX) where {TX <: XState} = _ctssyms(TX)

_ctssyms(TX::Type{<: XState}) = _syms(TX)[_findcts(TX)]

_ctssymstt(::Type{TX}) where {TX <: VState} = TX 

function _ctssymstt(::Type{TX}) where {TX <: PState}
   SYMS, TT = _symstt(TX)
   icts = _findcts(TX)
   CSYMS = SYMS[icts]
   CTT = Tuple( TT.types[i] for i in icts )
   return CSYMS, CTT
end

## ----- VState constructors

VState(t::NT) where {NT <: NamedTuple} = VState{NT}(t)

VState(; kwargs...) = VState(NamedTuple(kwargs))

@generated function VState(X::TX) where {TX <: PState}
   CSYMS = _ctssyms(TX) 
   quote 
      VState(select(_x(X), $CSYMS))
   end
end

      
"""
convert a PState to a corresponding VState 
(basically just remove the discrete variables)
"""
vstate_type(X::VState) = typeof(X)

@generated function vstate_type(X::TX)  where {TX <: PState}
   CSYMS = _ctssyms(TX) 
   quote
      typeof( VState( select(_x(X), $CSYMS) ) )
   end
end

vstate_type(::Type{TX}) where {TX <: VState} = TX 

# this shouldn't have to be generated, but that seemt to be the only 
# way I can avoid an allocation ? weird ? 
@generated function vstate_type(::Type{TX}) where {TX <: PState}
   CSYMS, TT = _ctssymstt(TX) 
   TTT = typeof( zero.(TT) )
   TV = VState{NamedTuple{CSYMS, TTT}}
   quote
      $TV
   end
end

# ------------ type promotion 

@generated function promote_rule(::Type{TDX}, ::Type{S}
                                ) where {TDX <: VState, S <: Number}
   SYMS, TT = _symstt(TDX)
   PTT = [ _mypromrl(S, TT.types[i]) for i = 1:length(SYMS) ]
   PTTstr = "Tuple{" * "$(tuple(PTT...))"[2:end-1] * "}"
   quote
      $(Meta.parse( "VState{NamedTuple{$(SYMS), $PTTstr}}" ))
   end
end


@generated function promote_rule(::Type{TDX1}, ::Type{TDX2}
                        ) where {TDX1 <: VState, TDX2 <: VState}
   SYMS1, TT1 = _symstt(TDX1)
   SYMS2, TT2 = _symstt(TDX2)
   SYMS = tuple(sort([union(SYMS1, SYMS2)...])...)
   PTT = [] 
   for sym in SYMS 
      if sym in SYMS1 && !(sym in SYMS2)
         push!(PTT, TT1.types[findfirst(isequal(sym), SYMS1)])
      elseif !(sym in SYMS1) && sym in SYMS2
         push!(PTT, TT2.types[findfirst(isequal(sym), SYMS2)])
      else 
         T1 = TT1.types[findfirst(isequal(sym), SYMS1)]
         T2 = TT2.types[findfirst(isequal(sym), SYMS2)]
         push!(PTT, _mypromrl(T1, T2))
      end
   end
   PTTstr = "Tuple{" * "$(tuple(PTT...))"[2:end-1] * "}"
   NTTex = Meta.parse("NamedTuple{$SYMS, $PTTstr}")
   quote
      return VState{$NTTex}
   end
end


# -------------------
      
# the next variant of vstate_type is used to potentially extend 
# from real states to complex vstates. 

_mypromrl(T::Type{<: Number}, S::Type{<: Number}) = 
      promote_type(T, S)
_mypromrl(T::Type{<: Number}, ::Type{<: SVector{N, P}}) where {N, P} = 
      SVector{N, promote_type(T, P)}
_mypromrl(::Type{<: SVector{N, P}}, T::Type{<: Number}) where {N, P} = 
      promote_type(T, P)
_mypromrl(T::Type{<: SVector{N, P1}}, ::Type{<: SVector{N, P2}}) where {N, P1, P2} = 
      SVector{N, promote_type(P1, P2)}

@generated function vstate_type(x::S, X::TX) where {S, TX <: XState}
   SYMS, TT = _symstt(TX)
   icts = _findcts(TX)
   CSYMS = SYMS[icts]
   CTT = [ _mypromrl(S, TT.types[i]) for i in icts ]
   CTTstr = "Tuple{" * "$(tuple(CTT...))"[2:end-1] * "}"
   quote
      $(Meta.parse( "VState{NamedTuple{$(CSYMS), $CTTstr}}" ))
   end
end

vstate_type(S::Type, X::XState) = vstate_type(zero(S), X)


## ---------- a simple copying setter-like functionality 

# this seems allocating, unclear why hence a generated implementation below 

@generated function setproperty(x::TX, s::Symbol, val) where {TX <: XState} 
   SYMS, TT = _symstt(TX)
   code = "$(nameof(TX))("
   for i = 1:length(SYMS)
      sym = SYMS[i]
      code *= "$sym = (:$sym == s ? val : x.$sym)::(typeof(x.$sym)), "
   end 
   code *= ")"
   quote
      $(Meta.parse(code))      
   end
end

@generated function setproperty(x::TX, ::Val{sym1}, val1) where {TX <: XState, sym1}
   SYMS, TT = _symstt(TX)
   code = "$(nameof(TX))("
   for (sym, T) in zip(SYMS, TT.types)
      if sym == sym1 
         code *= "$sym = val1, " 
      else 
         code *= "$sym = x.$sym, " 
      end 
   end
   code *= ")"
   quote
      $(Meta.parse(code))      
   end
end 


## ---------- explicit real/complex conversion 
# this feels a bit like a hack but might be unavoidable; 
# real, complex goes to _mod_real, _mod_complex, which is then applied 
# in only slightly non-standard fashion recursively to the states


for f in (:real, :imag, :complex, )
   fmod = Symbol("_mod_$f")
   eval(quote
      import Base: $f
      $fmod(x::Number) = $f(x)
      $fmod(x::StaticArrays.StaticArray) = $f.(x)
      function $f(X::TDX) where {TDX <: VState}
         SYMS = _syms(TDX)
         vals = ntuple(i -> $fmod(getproperty(X, SYMS[i])), length(SYMS))
         return TDX( NamedTuple{SYMS}(vals) )
      end
   end)
end

for f in (:real, :complex, )
   fmod = Symbol("_mod_$f")
   eval(quote
      $f(TDX::Type{<: VState}) = typeof( $f(zero(TDX)) )
      # import Base: $f
      # $fmod(x::Type{<: Number}) = $f(x)
      # $fmod(x::SVector{N, T}) where {N, T} = SVector{N, $f(T)}
      # function $f(TDX::Type{<: VState})
      #    SYMS, TT = _symstt(TDX)
      #    TT1 = ntuple(i -> $fmod(TT[i]), length(SYMS))
      #    vals = zero.(TT1)
      #    return typeof(TDX( NamedTuple{SYMS}(vals) )
      # end
   end)
end



for f in (:rand, :randn, :zero)
   fmod = Symbol("_mod_$f")
   eval( quote 
      import Base: $f 
      $fmod(T::Type) = $f(T) 
      $fmod(x::Union{Number, AbstractArray}) = $f(typeof(x))

      function $f(x::Union{TX, Type{TX}}) where {TX <: XState}
         SYMS, TT = _symstt(x)
         vals = ntuple(i -> $fmod(TT.types[i]), length(SYMS))
         return TX( NamedTuple{SYMS}( vals ) )
      end
   end )
end

# an extra for symbols, this is a bit questionable; why do we even need it?
_mod_zero(::Union{Symbol, Type{Symbol}}) = :O


## ----------- Some arithmetic operations 

# binary operations 

import Base: +, -

_prom_plus(TX1::Type{<: PState}, TX2::Type{<: PState}) = error("cannot add two PStates")
_prom_plus(TX1::Type{<: PState}, TX2::Type{<: VState}) = PState
_prom_plus(TX1::Type{<: VState}, TX2::Type{<: PState}) = PState
_prom_plus(TX1::Type{<: VState}, TX2::Type{<: VState}) = VState

function +(X1::TX1, X2::TX2) where {TX1 <: XState, TX2 <: XState}
   SYMS = _syms(X1, X2) 
   vals = ntuple( i -> begin 
            sym = SYMS[i]
            if haskey(_x(X1), sym) && haskey(_x(X2), sym)
               return getproperty(_x(X1), sym) + getproperty(_x(X2), sym)
            elseif haskey(_x(X1), sym)
               return getproperty(_x(X1), sym)
            else 
               return getproperty(_x(X2), sym)
            end
         end, length(SYMS))
   return _prom_plus(TX1, TX2)( NamedTuple{SYMS}(vals) )
end

_prom_minus(TX1::Type{<: PState}, TX2::Type{<: PState}) = VState
_prom_minus(TX1::Type{<: PState}, TX2::Type{<: VState}) = PState
_prom_minus(TX1::Type{<: VState}, TX2::Type{<: PState}) = error("Cannot subtrate a PState from a VState")
_prom_minus(TX1::Type{<: VState}, TX2::Type{<: VState}) = VState

function -(X1::TX1, X2::TX2) where {TX1 <: XState, TX2 <: XState}
   SYMS = _syms(X1, X2) 
   vals = ntuple( i -> begin 
            sym = SYMS[i]
            if haskey(_x(X1), sym) && haskey(_x(X2), sym)
               return getproperty(_x(X1), sym) - getproperty(_x(X2), sym)
            elseif haskey(_x(X1), sym)
               return getproperty(_x(X1), sym)
            else 
               return - getproperty(_x(X2), sym)
            end
         end, length(SYMS))
   return _prom_minus(TX1, TX2)( NamedTuple{SYMS}(vals) )
end


# multiplication with a scalar 
function *(X1::TX, a::Number) where {TX <: XState}
   SYMS = _syms(TX)
   vals = ntuple( i -> *( getproperty(_x(X1), SYMS[i]), a ), length(SYMS) )
   return TX( NamedTuple{SYMS}(vals) )
end

*(a::Number, X1::XState) = *(X1, a)

*(aa::SVector{N, <: Number}, X1::XState) where {N} = aa .* Ref(X1)
promote_rule(::Type{SVector{N, T}}, ::Type{TX}) where {N, T <: Number, TX <: XState} = 
      SVector{N, promote_type(T, TX)}

# unary 
import Base: - 

for f in (:-, )
   eval(quote
      function $f(X::TX) where {TX <: XState}
         SYMS = _syms(TX)
         vals = ntuple( i -> $f( getproperty(_x(X), SYMS[i]) ), length(SYMS) )
         return TX( NamedTuple{SYMS}(vals) )
      end
   end)
end


# reduction to scalar 

import LinearAlgebra: dot 
import Base: isapprox



for (f, g) in ( (:dot, :sum), (:isapprox, :all) )  # (:contract, :sum), 
   eval( quote 
      function $f(X1::TX1, X2::TX2) where {TX1 <: XState, TX2 <: XState}
         SYMS = _syms(TX1)
         @assert SYMS == _syms(TX2)
         return $g( $f( getproperty(_x(X1), sym), 
                         getproperty(_x(X2), sym) )   for sym in SYMS)
      end
   end )
end

@generated function contract(X1::TX1, X2::TX2) where {TX1 <: XState, TX2 <: XState}
   # this line is important - it means that missing symbols are interpreted as zero
   SYMS = intersect(_syms(TX1), _syms(TX2))
   code = "contract(X1.$(SYMS[1]), X2.$(SYMS[1]))"
   for sym in SYMS[2:end]
      code *= " + contract(X1.$sym, X2.$sym)"
   end
   return quote 
      $(Meta.parse(code))
   end
end

contract(X1::Number, X2::XState) = X1 * X2 
contract(X2::XState, X1::Number) = X1 * X2 
contract(x1::Number, x2::Number) = x1 * x2 
contract(x1::SVector, x2::SVector) = sum(x1 .* x2)

import LinearAlgebra: norm 

for (f, g) in ((:norm, :norm), (:sumsq, :sum), (:normsq, :sum) )
   eval( quote 
      function $f(X::TX) where {TX <: XState}
         SYMS = _syms(TX)
         vals = ntuple( i -> $f( getproperty(_x(X), SYMS[i]) ),  length(SYMS))
         return $g(vals)
      end

      # this is a corner case which is a bit unclear but for now it seems 
      # convenient to have this 
      $f(::VState{@NamedTuple{}}) = zero(Bool)
      $f(::PState{@NamedTuple{}}) = zero(Bool)
   end )
end


sumsq(x::SVector) = sum(x.^2)
normsq(x::SVector) = dot(x, x)




## --------  not clear where needed; or deleted functionality 

# function promote_leaf_eltypes(X::TX) where {TX <: XState} 
#    SYMS = _syms(TX)
#    promote_type( ntuple(i -> promote_leaf_eltypes(getproperty(_x(X), SYMS[i])), length(SYMS))... )
# end



# ---------------- AD code 

# TODO: check whether this is still needed 
# this function makes sure that gradients w.r.t. a PState become a VState 
# The original iple,entation is a bit more complicated, but I don't understand 
# it anymore. I'll need to revisit this a bit. 
# function rrule(::typeof(getproperty), X::XState, sym::Symbol) 
#    val = getproperty(X, sym)
#    return val, w -> ( NoTangent(), 
#                       vstate_type(w[1], X)( NamedTuple{(sym,)}((w,)) ), 
#                       NoTangent() )
# end


import ChainRulesCore
import ChainRulesCore: rrule, NoTangent 

function rrule(::typeof(Base.getproperty), X::XState, sym::Symbol)
   val = getproperty(X, sym)
   return val, Δ -> (NoTangent(), VState(NamedTuple{(sym,)}((Δ,))), NoTangent())
end
