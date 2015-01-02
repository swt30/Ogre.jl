module Common
# Exported variables
export DATADIR
# Exported types
export Callable, Equation, EquationSet, ValueSet
# Exported functions
export callfunc, zero

abstract Callable
abstract Equation <: Callable

type EquationSet{T<:Equation} <: Callable
    equations::Vector{T}
end

n_eqs(es::EquationSet) = length(es.equations)

type ValueSet{T<:Real}
    m::T
    r::T
    P::T
end

import Base.zero
zero(::Type{ValueSet}) = ValueSet(0, 0, 0)

function callfunc(eq::Equation, x::Real)
    eq.equation(x)
end

function callfunc(eq::Equation, vs::ValueSet)
    eq.equation(vs)
end

function callfunc(es::EquationSet, vs::ValueSet)
    [callfunc(equation, vs) for equation in es.equations]
end

function callfunc{T<:Real}(cl::Callable, t::T, y::Vector{T})
    callfunc(cl, ValueSet(t, y...))
end

const DATADIR = Pkg.dir("Ogre", "data")

end
