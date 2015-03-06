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

import Base.length
length(es::EquationSet) = length(es.equations)

# Type for passing around physical parameters
#------------------------------------------------------------------------------

abstract ModelComplexity
immutable NoTemp <: ModelComplexity; end
immutable WithTemp <: ModelComplexity; end

const notemp = NoTemp()
const withtemp = WithTemp()

abstract ValueSet

@doc """Holds physical values of mass, radius, pressure, and temperature.""" ->
immutable PhysicalValues{R<:Real} <: ValueSet
    m::R
    r::R
    P::R
    T::R
end
function PhysicalValues(m::Real, r::Real, P::Real, T::Real)
    PhysicalValues(promote(m, r, P, T)...)
end

@doc """Holds physical values of mass, radius, and pressure""" ->
immutable MassRadiusPressure{R<:Real} <: ValueSet
    m::R
    r::R
    P::R
end
function MassRadiusPressure(m::Real, r::Real, P::Real)
    MassRadiusPressure(promote(m, r, P)...)
end

# Constructors
ValueSet(m, r, P) = ValueSet(notemp, m, r, P)
ValueSet(m, r, P, T) = ValueSet(withtemp, m, r, P, T)
ValueSet(::NoTemp, args...) = MassRadiusPressure(args...)
ValueSet(::WithTemp, args...) = PhysicalValues(args...)

# Properties
@doc "Get the dependent physical values (radius, pressure, [temperature])" ->
nonmass(pv::PhysicalValues) = [pv.r, pv.P, pv.T]
nonmass(mrp::MassRadiusPressure) = [mrp.r, mrp.P]
@doc "Get the independent physical coordinate (mass)" ->
mass(vs::ValueSet) = vs.m
radius(vs::ValueSet) = vs.r
pressure(vs::ValueSet) = vs.P
temperature(vs::PhysicalValues) = vs.T

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
zero(::Type{ValueSet}, ::NoTemp) = zero(MassRadiusPressure)
zero(::Type{ValueSet}, ::WithTemp) = zero(PhysicalValues)
call(eq::Equation, x::Real) = eq.equation(x)
call(eq::Equation, vs::ValueSet) = eq.equation(vs)
call(es::EquationSet, vs::ValueSet) =  map(eq -> eq(vs), es.equations)
call{T<:Real}(cl::Callable, t::T, y::Vector{T}) = cl(ValueSet(t, y...))

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
