# Types and functions for passing around sets of physical values


# Types

"Holds a set of planetary physical values (mass, radius, etc)"
abstract ValueSet{mc<:ModelComplexity}

"Holds physical values of mass, radius, pressure, and temperature."
type PhysicalValues{R<:Real} <: ValueSet{WithTemp}
    m::R
    r::R
    P::R
    T::R
end
PhysicalValues(m, r, P, T) = PhysicalValues(promote(m, r, P, T)...)
ValueSet(m, r, P, T) = PhysicalValues(m, r, P, T)

"Holds physical values of mass, radius, and pressure"
type MassRadiusPressure{R<:Real} <: ValueSet{NoTemp}
    m::R
    r::R
    P::R
end
MassRadiusPressure(m, r, P) = MassRadiusPressure(promote(m, r, P)...)
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
    m = m / M_earth
    r = r / R_earth
    println("$m M⊕, $r R⊕, $P Pa, $T K")
end
function Base.show(io::IO, mrp::MassRadiusPressure)
    m, r, P = mrp.m, mrp.r, mrp.P
    m = m / M_earth
    r = r / R_earth
    println("$m M⊕, $r R⊕, $P Pa")
end

# Get individual values from the ValueSets
mass(vs::ValueSet) = vs.m
radius(vs::ValueSet) = vs.r
pressure(vs::ValueSet) = vs.P
temperature(vs::PhysicalValues) = vs.T
function gravity(vs::ValueSet)
    m = mass(vs)
    r = radius(vs)
    G * m ./ (r.^2)
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
