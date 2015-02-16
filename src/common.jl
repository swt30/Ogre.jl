import Base.zero, Base.call

#= Callable and equation types =#
abstract Callable
abstract Equation <: Callable

@doc """
    A set of equations which can be evaluated all at once.

    `equations`: Vector of `Equation`
    """ ->
immutable EquationSet{T<:Equation} <: Callable
    equations::Vector{T}
end

n_eqs(es::EquationSet) = length(es.equations)

@doc """Holds physical values of mass, radius, and pressure.""" ->
immutable ValueSet{T<:Real}
    m::T
    r::T
    P::T
end

@doc "Get the dependent physical values (radius, pressure)" ->
physical_values{T<:Real}(vs::ValueSet{T}) = [vs.r, vs.P]::Vector{T}
@doc "Get the independent physical coordinate (mass)" ->
mass_coordinate(vs::ValueSet) = vs.m::Real

zero(::Type{ValueSet}) = ValueSet(0, 0, 0)
call(eq::Equation, x::Real) = eq.equation(x)
call(eq::Equation, vs::ValueSet) = eq.equation(vs)
call(es::EquationSet, vs::ValueSet) =  map(eq -> eq(vs), es.equations)
call{T<:Real}(cl::Callable, t::T, y::Vector{T}) = cl(ValueSet(t, y...))

@doc """Location of the data files in this package""" ->
const DATADIR = Pkg.dir("Ogre", "data")

@doc """Copy a type, modifying certain fields""" ->
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
# TODO: Remove uses of cpmod once default constructors for immutables are
# changed
