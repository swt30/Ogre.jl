# COMMON.JL
# Common functionality across all of Ogre.jl

# Callable and equation types
#-------------------------------------------------------------------------------

"Can be called using the standard call syntax"
abstract Callable
"Represents a physical equation"
abstract Equation <: Callable
"Represents a set of equations which can be evaluated all at once"
immutable EquationSet <: Callable
    equations::Vector{Equation}
end

Base.length(es::EquationSet) = length(es.equations)
Base.getindex(es::EquationSet, i...) = es.equations[i...]

# Types for passing around physical parameters
#------------------------------------------------------------------------------

"Define the complexity of a physical system"
abstract ModelComplexity
"This system excludes temperature details"
immutable NoTemp <: ModelComplexity; end
"This system explicitly includes temperature details"
immutable WithTemp <: ModelComplexity; end
"This system explicitly includes both temperature and pressure details"
immutable WithTempPressure <: ModelComplexity; end

"Number of physical variables (mass, radius, pressure...) used"
function nvars end
nvars(::Type{NoTemp}) = 3
nvars(::Type{WithTemp}) = 4

"Number of independent physical variables (radius, pressure...) used"
ndeps{mc<:ModelComplexity}(::Type{mc}) = nvars(mc) - 1

"Holds a set of planetary physical values (mass, radius, etc)"
abstract ValueSet{mc<:ModelComplexity}

"Holds physical values of mass, radius, pressure, and temperature."
immutable PhysicalValues{R<:Real} <: ValueSet{WithTemp}
    m::R
    r::R
    P::R
    T::R
end
PhysicalValues(m, r, P, T) = PhysicalValues(promote(m, r, P, T)...)
ValueSet(m, r, P, T) = PhysicalValues(m, r, P, T)

function Base.show(io::IO, pv::PhysicalValues)
    m, r, P, T = pv.m, pv.r, pv.P, pv.T
    m = m / M_earth
    r = r / R_earth
    println("$m M⊕, $r R⊕, $P Pa, $T K")
end

"Holds physical values of mass, radius, and pressure"
immutable MassRadiusPressure{R<:Real} <: ValueSet{NoTemp}
    m::R
    r::R
    P::R
end
MassRadiusPressure(m, r, P) = MassRadiusPressure(promote(m, r, P)...)
ValueSet(m, r, P) = MassRadiusPressure(m, r, P)

function Base.show(io::IO, mrp::MassRadiusPressure)
    m, r, P = pv.m, pv.r, pv.P
    m = m / M_earth
    r = r / R_earth
    println("$m M⊕, $r R⊕, $P Pa")
end

# Properties of ValueSets
mass(vs::ValueSet) = vs.m
radius(vs::ValueSet) = vs.r
pressure(vs::ValueSet) = vs.P
temperature(vs::PhysicalValues) = vs.T
function gravity(vs::ValueSet)
    m = mass(vs)
    r = radius(vs)
    @show G * m ./ (r.^2)
end
"Get dependent physical values (radius, pressure, [temperature])"
function nonmass end
nonmass(pv::PhysicalValues) = [radius(pv), pressure(pv), temperature(pv)]
nonmass(mrp::MassRadiusPressure) = [radius(mrp), pressure(mrp)]

"Is a given `ValueSet` physical (all positive)?"
function isphysical end
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
"Location of the data files in this package"
const DATADIR = Pkg.dir("Ogre", "data")

"Copy a type, modifying specified fields"
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

"Create an interpolating function on a linear grid"
function lininterp(xs::Vector, ys::Vector)
    spline = Spline1D(xs, ys, k=2)

    interp_func(x::Real) = evaluate(spline, Float64(x))
    interp_func(x::AbstractVector) = evaluate(spline, Vector{Float64}(x))

    interp_func
end
function lininterp{S, T}(xs::AbstractVector{S}, ys::AbstractVector{T})
    lininterp(Vector{S}(xs), Vector{T}(ys))
end

"Create an interpolating function using a log-spaced coordinate grid."
function loginterp(xs::Vector, ys::Vector)
    # first transform the grid to be linear
    logxs = log10(xs)
    # then do the interpolation as if it were linear
    lin_interp_func = lininterp(logxs, ys)

    interp_func(x) = lin_interp_func(log10(x))
end
function loginterp{S, T}(xs::AbstractVector{S}, ys::AbstractVector{T})
    loginterp(Vector{S}(xs), Vector{T}(ys))
end

"Create an interpolating function from a 2D linear grid"
function lininterp(xs::Vector, ys::Vector, zs::Matrix)
    hasnan(zs) ? error("2D data contains NaNs") :
    spline = Spline2D(xs, ys, zs, kx=1, ky=1)

    interp_func(x::Real, y::Real) = evaluate(spline, Float64(x), Float64(y))
    function interp_func(x::AbstractVector, y::AbstractVector)
        evaluate(spline, Vector{Float64}(x), Vector{Float64}(y))
    end

    interp_func
end
function lininterp{S, T}(xs::AbstractVector{S}, ys::AbstractVector{T}, zs; kwargs...)
    lininterp(Vector{S}(xs), Vector{T}(ys), zs; kwargs...)
end

"Create an interpolating function from a 2D log-log grid"
function loginterp(xs::Vector, ys::Vector, zs::Matrix; kwargs...)
    logxs = log10(xs)
    logys = log10(ys)

    lin_interp_func = lininterp(logxs, logys, zs; kwargs...)
    interp_func(x, y) = lin_interp_func(log10(x), log10(y))
end

"Create an interpolating function from a 2D log-linear grid"
function semiloginterpx(xs::Vector, ys::Vector, zs::Matrix; kwargs...)
    logxs = log10(xs)

    lin_interp_func = lininterp(logxs, ys, zs; kwargs...)
    interp_func(x, y) = lin_interp_func(log10(x), y)
end

# Miscellaneous utility funcs
#-------------------------------------------------------------------------------

"Map a function across rows of a 2D array"
maprows(f, m::Matrix) = mapslices(f, m, 2)
"Does this contain any NaN values?"
hasnan(x) = any(isnan(x))
notnan(x) = !isnan(x)

"Wait for ENTER to be pressed before continuing"
function wait_for_enter()
    println("Press ENTER to continue...")
    readline(STDIN)
end