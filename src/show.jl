using Unitful: Quantity, @u_str, unit, ustrip

const _showdigits = Ref{Int64}(2)

"""
change how many digits are printed in the `PState` or `VState` output 
"""
function showdigits!(n::Integer)
   _showdigits[] = n 
end 

_2str(x) = string(x)

_2str(x::AbstractFloat) = "$(round(x, digits=_showdigits[]))"
_2str(x::Complex) = "$(round(x, digits=_showdigits[]))"
_2str(x::Quantity) = _2str(ustrip(x)) * " $(unit(x))"

_2str(x::SVector{N, <: AbstractFloat}) where {N} = string(round.(x, digits=_showdigits[]))
_2str(x::SVector{N, <: Complex}) where {N} = string(round.(x, digits=_showdigits[]))[11:end]
_2str(x::SVector{N, <: Quantity}) where {N} = _2str(ustrip.(x)) * " $(unit(x[1]))"

# MATHEMATICAL LEFT WHITE SQUARE BRACKET (U+27E6, Ps): ⟦
# MATHEMATICAL RIGHT WHITE SQUARE BRACKET (U+27E7, Pe): ⟧
# MATHEMATICAL LEFT WHITE TORTOISE SHELL BRACKET (U+27EC, Ps): ⟬
# MATHEMATICAL RIGHT WHITE TORTOISE SHELL BRACKET (U+27ED, Pe): ⟭
# LEFT DOUBLE PARENTHESIS (U+2E28, Ps): ⸨
# RIGHT DOUBLE PARENTHESIS (U+2E29, Pe): ⸩
# LEFT DOUBLE ANGLE BRACKET (U+300A, Ps): 《
# RIGHT DOUBLE ANGLE BRACKET (U+300B, Pe): 》
# LEFT WHITE LENTICULAR BRACKET (U+3016, Ps): 〖
# RIGHT WHITE LENTICULAR BRACKET (U+3017, Pe): 〗
# LEFT WHITE SQUARE BRACKET (U+301A, Ps): 〚
# RIGHT WHITE SQUARE BRACKET (U+301B, Pe): 〛
# FULLWIDTH LEFT WHITE PARENTHESIS (U+FF5F, Ps): ｟
# FULLWIDTH RIGHT WHITE PARENTHESIS (U+FF60, Pe): ｠

_leftbra( X::PState) = "〖"
_rightbra(X::PState) = "〗"
_leftbra( X::VState) = "｟"
_rightbra(X::VState) = "｠"

function Base.show(io::IO, X::XState)
   _str(sym) = "$(sym)=$(_2str(getproperty(_x(X), sym)))"
   strs = [ _str(sym) for sym in keys(_x(X)) ]
   str = prod(strs[i] * ", " for i = 1:length(strs)-1; init="") * strs[end] 
   # str = prod( "$(sym):$(_2str(getproperty(_x(X), sym))), " 
   #             for sym in keys(_x(X)) )
   print(io,  _leftbra(X) * str * _rightbra(X))
end


function Base.show(io::IO, sys::Union{AosSystem, SoaSystem})
   println(io, typeof(sys).name.name, ": len = $(length(sys)), ") 
   print(io, "  ", cell(sys))
   for i = 1:min(4, length(sys))
      println(io, "  ", sys[i])
   end
   if length(sys) > 4
      println(io, "  ... $(length(sys)-4) particles not shown ...")
   end
end

Base.show(io::IO, ::MIME"text/plain", sys::Union{AosSystem, SoaSystem}) = 
      show(io, sys)
