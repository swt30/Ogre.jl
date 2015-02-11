module Common
# Exported variables
export DATADIR
# Exported types
export Callable, Equation, EquationSet, ValueSet
# Exported functions
export callfunc, zero, cpmod

#= Callable and equation types =#

abstract Callable
abstract Equation <: Callable

immutable EquationSet{T<:Equation} <: Callable
    equations::Vector{T}
end

n_eqs(es::EquationSet) = length(es.equations)

immutable ValueSet{T<:Real}
    m::T
    r::T
    P::T
end

Base.zero(::Type{ValueSet}) = ValueSet(0, 0, 0)

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

#= Location of the data files =#

const DATADIR = Pkg.dir("Ogre", "data")

#= A useful function to copy and modify an immutable type =#

function cpmod{T}(pp::T, di)
    di = !isa(di, Associative) ? Dict(di) : di
    ns = names(pp)
    args = Array(Any, length(ns))
    for (i,n) in enumerate(ns)
        args[i] = get(di, n, getfield(pp, n))
    end
    T(args...)
end
cpmod{T}(pp::T; kws...) = cpmod(pp, kws)

end
