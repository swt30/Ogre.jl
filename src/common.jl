module common
export Callable, Equation, EquationSet, ValueSet, call, datadir

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

function call(eq::Equation, x::Float64)
    eq.equation(x)
end

function call(eq::Equation, vs::ValueSet)
    eq.equation(vs)
end

function call(es::EquationSet, vs::ValueSet)
    [call(equation, vs) for equation in es.equations]
end

function call(cl::Callable, t::Float64, y::Vector{Float64})
    call(cl, ValueSet(t, y...))
end

datadir() = Pkg.dir("ogre", "data")

end
