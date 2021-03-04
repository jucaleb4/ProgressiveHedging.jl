"""
Abstract type for ProgressiveHedging penalty parameter.

Concrete sybtypes determine how this penalty is used within the PH algorithm.
"""
abstract type AbstractPenaltyParameter end

## Scalar Penalty Implementation ## 

struct ScalarPenaltyParameter <: AbstractPenaltyParameter
    value::Float64
end

## Proportional Penalty Implementation ##
struct ProportionalPenaltyParameter <: AbstractPenaltyParameter
    constant::Float64
    coefficient::Dict{XhatID,Float64}
end

function ProportionalPenaltyParameter(constant::Float64)
    return ProportionalPenaltyParameter(constant, Dict{XhatID,Float64}())
end

# NOTE: temporary fcns for coefficients of expressions
# Future JuMP versions will implement this
coefficient(a::JuMP.GenericAffExpr{C,V}, v::V) where {C,V} = get(a.terms, v, zero(C))
coefficient(a::JuMP.GenericAffExpr{C,V}, v1::V, v2::V) where {C,V} = zero(C)
function coefficient(q::JuMP.GenericQuadExpr{C,V}, v1::V, v2::V) where {C,V}
    return get(q.terms, UnorderedPair(v1,v2), zero(C))
end
coefficient(q::JuMP.GenericQuadExpr{C,V}, v::V) where {C,V} = coefficient(q.aff, v)

function penalty_value(r::AbstractPenaltyParameter, var)
    throw(UnimplementedError("penalty_value is unimplemented for `r` of type $(typeof(r)) and `var` of type $(typeof(var))."))
end

function penalty_value(r::ScalarPenaltyParameter, 
                        obj::JuMP.GenericQuadExpr,
                        var::JuMP.VariableRef
                        )::Float64
    return r.value
end

function penalty_value(r::ProportionalPenaltyParameter, 
                        obj::JuMP.GenericQuadExpr,
                        var::JuMP.VariableRef
                        )::Float64
    coeff = coefficient(obj, var)
    return r.constant * coeff
end

function penalty_value(r::ProportionalPenaltyParameter,
                       xhid::XhatID
                       )::Float64
    return r.constant * r.coefficient[xhid]
end