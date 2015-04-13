# COMMON.JL
# Common functionality across all of Ogre

# Callable and equation types

abstract Callable
abstract Equation <: Callable

@doc """A set of equations which can be evaluated all at once.

    `equations`: Vector of `Equation` """ ->
immutable EquationSet <: Callable
    equations::Vector{Equation}
end

import Base.length
length(es::EquationSet) = length(es.equations)

# Type for passing around physical parameters
#------------------------------------------------------------------------------

abstract ModelComplexity
immutable NoTemp <: ModelComplexity; end
immutable WithTemp <: ModelComplexity; end

@doc "Number of physical variables used" ->
nvars(::Type{NoTemp}) = 3
nvars(::Type{WithTemp}) = 4
@doc "Number of independent physical variables used" ->
ndeps{mc<:ModelComplexity}(::Type{mc}) = nvars(mc) - 1

abstract ValueSet{mc<:ModelComplexity}

@doc "Holds physical values of mass, radius, pressure, and temperature." ->
immutable PhysicalValues{R<:Real} <: ValueSet{WithTemp}
    m::R
    r::R
    P::R
    T::R
end
function PhysicalValues(m::Real, r::Real, P::Real, T::Real)
    PhysicalValues(promote(m, r, P, T)...)
end
ValueSet(m, r, P, T) = PhysicalValues(m, r, P, T)


@doc "Holds physical values of mass, radius, and pressure" ->
immutable MassRadiusPressure{R<:Real} <: ValueSet{NoTemp}
    m::R
    r::R
    P::R
end
function MassRadiusPressure(m::Real, r::Real, P::Real)
    MassRadiusPressure(promote(m, r, P)...)
end
ValueSet(m, r, P) = MassRadiusPressure(m, r, P)

# Properties
@doc "Get independent physical coordinate (mass)" ->
mass(vs::ValueSet) = vs.m
@doc "Get radius coordinate" ->
radius(vs::ValueSet) = vs.r
@doc "Get the pressure" ->
pressure(vs::ValueSet) = vs.P
@doc "Get the temperature" ->
temperature(vs::PhysicalValues) = vs.T
@doc "Get dependent physical values (radius, pressure, [temperature])" ->
nonmass(pv::PhysicalValues) = [radius(pv), pressure(pv), temperature(pv)]
nonmass(mrp::MassRadiusPressure) = [radius(mrp), pressure(mrp)]

@doc "Is a given `ValueSet` physical (all positive?)"
function isphysical(vs::PhysicalValues)
    mass(vs) > 0 && radius(vs) > 0 && pressure(vs) > 0 && temperature(vs) > 0
end
function isphysical(vs::MassRadiusPressure)
    mass(vs) > 0 && radius(vs) > 0 && pressure(vs) > 0
end

# Interaction
import Base: zero, call
zero(::Type{MassRadiusPressure}) = MassRadiusPressure(0, 0, 0)
zero(::Type{PhysicalValues}) = PhysicalValues(0, 0, 0, 0)
zero(::Type{ValueSet{NoTemp}}) = zero(MassRadiusPressure)
zero(::Type{ValueSet{WithTemp}}) = zero(PhysicalValues)
call(eq::Equation, x::Real) = eq.equation(x)
call(eq::Equation, vs::ValueSet) = eq.equation(vs)
call(es::EquationSet, vs::ValueSet) =  map(eq -> eq(vs), es.equations)
call{T<:Real}(cl::Callable, t::T, y::Vector{T}) = cl(ValueSet(t, y...))

# Useful common definitions
@doc "Location of the data files in this package" ->
const DATADIR = Pkg.dir("Ogre", "data")

@doc "Copy a type, modifying certain fields" ->
function cpmod{T}(pp::T, di)
    di = !isa(di, Associative) ? Dict(di) : di
    ns = fieldnames(pp)
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
