module common
export Callable, Equation, EquationSet, ValueSet, callfunc, datadir, zero
import Base.zero

abstract Callable
abstract Equation <: Callable

type EquationSet{T<:Equation} <: Callable
    equations::Vector{T}
end

n_eqs(es::EquationSet) = length(es.equations)

type ValueSet
    m::Float64
    r::Float64
    P::Float64
end

zero(::Type{ValueSet}) = ValueSet(0, 0, 0)

function callfunc(eq::Equation, x::Float64)
    eq.equation(x)
end

function callfunc(eq::Equation, vs::ValueSet)
    eq.equation(vs)
end

function callfunc(es::EquationSet, vs::ValueSet)
    [callfunc(equation, vs) for equation in es.equations]
end

function callfunc(cl::Callable, t::Float64, y::Vector{Float64})
    callfunc(cl, ValueSet(t, y...))
end

datadir() = Pkg.dir("ogre", "data")

end
