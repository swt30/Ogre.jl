# COMMON.JL
# Common functionality across all of Ogre

# Callable and equation types

abstract Callable
abstract Equation <: Callable

@doc """
    A set of equations which can be evaluated all at once.

    `equations`: Vector of `Equation`
    """ ->
immutable EquationSet{T<:Equation} <: Callable
    equations::Vector{T}
end

Base.length(es::EquationSet) = length(es.equations)

# Type for passing around physical parameters
#------------------------------------------------------------------------------

abstract ValueSet

@doc """Holds physical values of mass, radius, pressure, and temperature.""" ->
immutable PhysicalValues{E<:Real} <: ValueSet
    m::E
    r::E
    P::E
    T::E
end
function PhysicalValues(m::Real, r::Real, P::Real, T::Real)
    PhysicalValues(promote(m, r, P, T)...)
end

@doc """Holds physical values of mass, radius, and pressure""" ->
immutable MassRadiusPressure{T<:Real} <: ValueSet
    m::T
    r::T
    P::T
end
function MassRadiusPressure(m::Real, r::Real, P::Real)
    MassRadiusPressure(promote(m, r, P)...)
end

ValueSet(m::Real, r::Real, P::Real) = MassRadiusPressure(m, r, P)
ValueSet(m::Real, r::Real, P::Real, T::Real) = PhysicalValues(m, r, P, T)

@doc "Get the dependent physical values (radius, pressure, [temperature])" ->
nonmass(pv::PhysicalValues) = [pv.r, pv.P, pv.T]
nonmass(mrp::MassRadiusPressure) = [mrp.r, mrp.P]
@doc "Get the independent physical coordinate (mass)" ->
mass(vs::ValueSet) = vs.m

Base.zero(::Type{PhysicalValues}) = PhysicalValues(0, 0, 0, 0)
Base.zero(::Type{MassRadiusPressure}) = MassRadiusPressure(0, 0, 0)
Base.call(eq::Equation, x::Real) = eq.equation(x)
Base.call(eq::Equation, vs::ValueSet) = eq.equation(vs)
Base.call(es::EquationSet, vs::ValueSet) =  map(eq -> eq(vs), es.equations)
Base.call{T<:Real}(cl::Callable, t::T, y::Vector{T}) = cl(ValueSet(t, y...))

# Useful common definitions
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

# Physical constants

include("constants.jl")
