# Types and functions for passing around sets of physical values


# Types

"Holds a set of planetary physical values (mass, radius, etc)"
abstract type ValueSet{mc<:ModelComplexity} end

"Holds physical values of mass, radius, pressure, and temperature."
mutable struct PhysicalValues{Ma<:Mass, Ra<:Radius, Pr<:Pressure, Te<:Temperature} <: ValueSet{WithTemp}
    m::Ma
    r::Ra
    P::Pr
    T::Te
end
ValueSet(m, r, P, T) = PhysicalValues(m, r, P, T)

"Holds physical values of mass, radius, and pressure"
mutable struct MassRadiusPressure{Ma<:Mass,Ra<:Radius,Pr<:Pressure} <: ValueSet{NoTemp}
    m::Ma
    r::Ra
    P::Pr
end
ValueSet(m, r, P) = MassRadiusPressure(m, r, P)


# Query and work with sets of physical values

"Number of physical variables (mass, radius, pressure...) used"
function nvars end
nvars(::Type{NoTemp}) = 3  # M, R, P
nvars(::Type{WithTemp}) = 4  # M, R, P, T

"Number of independent physical variables (radius, pressure...) used"
ndeps{mc<:ModelComplexity}(::Type{mc}) = nvars(mc) - 1

# Print them in terms of Earth masses and radii
function Base.show(io::IO, pv::PhysicalValues)
    m, r, P, T = pv.m, pv.r, pv.P, pv.T
    show(io, "$(m/M_earth) M⊕, $(r/R_earth) R⊕, $(P/Pa) Pa, $(T/K) K")
end
function Base.show(io::IO, mrp::MassRadiusPressure)
    m, r, P = mrp.m, mrp.r, mrp.P
    show(io, "$(m/M_earth) M⊕, $(r/R_earth) R⊕, $(P/Pa) Pa")
end

# Get individual values from the ValueSets
mass(vs::ValueSet) = vs.m::Mass
radius(vs::ValueSet) = vs.r::Radius
pressure(vs::ValueSet) = vs.P::Pressure
temperature(vs::PhysicalValues) = vs.T::Temperature
function gravity(vs::ValueSet)
    m = mass(vs)
    r = radius(vs)
    (G * m ./ (r.^2))::Gravity
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

# Zero values
Base.zero(::Type{MassRadiusPressure}) = MassRadiusPressure(0, 0, 0)
Base.zero(::Type{PhysicalValues}) = PhysicalValues(0, 0, 0, 0)
Base.zero(::Type{ValueSet{NoTemp}}) = zero(MassRadiusPressure)
Base.zero(::Type{ValueSet{WithTemp}}) = zero(PhysicalValues)

# Set individual values from the ValueSets
setmass!(vs::ValueSet, m) = (vs.m = m)
setradius!(vs::ValueSet, r) = (vs.r = r)
setpressure!(vs::ValueSet, P) = (vs.P = P)
settemperature!(vs::PhysicalValues, T) = (vs.T = T)
