## ----- Implementation of a Position State, as a basic example 

PositionState{T} = State{NamedTuple{(:rr,), Tuple{SVector{3, T}}}}

PositionState(r::AbstractVector{T}) where {T <: AbstractFloat} = 
      (@assert length(r) == 3; PositionState{T}(; rr = SVector{3, T}(r)))

promote_rule(::Union{Type{S}, Type{PositionState{S}}}, 
             ::Type{PositionState{T}}) where {S, T} = 
      PositionState{promote_type(S, T)}

# some special functionality for PositionState, mostly needed for testing  
*(A::AbstractMatrix, X::TX) where {TX <: PositionState} = TX( (rr = A * X.rr,) )
+(X::TX, u::StaticVector{3}) where {TX <: PositionState} = TX( (rr = X.rr + u,) )


