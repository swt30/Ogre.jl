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
function PhysicalValues(m, r, P, T)
    PhysicalValues(promote(m, r, P, T)...)
end
ValueSet(m, r, P, T) = PhysicalValues(m, r, P, T)


@doc "Holds physical values of mass, radius, and pressure" ->
immutable MassRadiusPressure{R<:Real} <: ValueSet{NoTemp}
    m::R
    r::R
    P::R
end
function MassRadiusPressure(m, r, P)
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
call(cl::Callable, t::Real, y::Vector) = cl(ValueSet(t, y...))

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

# Interpolation
#-------------------------------------------------------------------------------

@doc "Create an interpolating function in linear space" ->
function lininterp(xs, ys)
    spline = Spline1D(xs, ys, k=2)
    interp_func(x::Real) = evaluate(spline, Float64(x))
    interp_func(x::Vector) = evaluate(spline, Vector{Float64}(x))

    interp_func
end

@doc "Create an interpolating function using a log-spaced coordinate grid." ->
function loginterp(xs, ys)
    # first transform the grid to be linear
    logxs = log10(xs)
    # then do the interpolation as if it were linear
    lin_interp_func = lininterp(logxs, ys)

    interp_func(x) = lin_interp_func(log10(x))
end

@doc "Create a linear interpolating function from a 2D grid" ->
function lininterp(xs, ys, zs::Matrix; suppress_warnings=false)
    if hasnan(zs)
        if !suppress_warnings
            warn("2D data contains NaNs: setting to sentinel value of -1e99")
        end
        zs[isnan(zs)] = -1e99
    end
    spline = Spline2D(xs, ys, zs, kx=1, ky=1)
    interp_func(x::Real, y::Real) = evaluate(spline, Float64(x), Float64(y))
    interp_func(x::Vector, y::Vector) = evaluate(spline, x, y)

    interp_func
end

@doc "Create a log-spaced interpolating function from a 2D grid" ->
function loginterp(xs, ys, zs::Matrix; kwargs...)
    logxs = log10(xs)
    logys = log10(ys)

    lin_interp_func = lininterp(logxs, logys, zs; kwargs...)
    interp_func(x, y) = lin_interp_func(log10(x), log10(y))
end
