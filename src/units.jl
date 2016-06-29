# Definitions of physical unit types for dispatch. At the moment this is a
# "dummy" set of definitions, but this may be swapped out for actual units if I
# have time.

# Put everything in SI units
const m = 1
const kg = 1
const s = 1
const K = 1
# Derived units
const W = kg*m^2/s^3
const Pa = kg/m/s^2
const J = kg*m^2/s^2
# Multiples of base units
const cm = 1e-2m
const g = 1e-3kg
const bar = 10^5 * Pa
const GPa = 10^9 * Pa

#= definitions unused for now
using SIUnits
using SIUnits.ShortUnits
import SIUnits.unit, SIUnits.SIQuantity

# typealias Unit{T} SIQuantity{T,...}
typealias Mass{T} SIQuantity{T,0,1,0,0,0,0,0,0,0}
typealias Distance{T} SIQuantity{T,1,0,0,0,0,0,0,0,0}
typealias Radius Distance
typealias Temperature{T} SIQuantity{T,0,0,0,0,1,0,0,0,0}
typealias MassLuminosity{T} SIQuantity{T,2,0,-3,0,0,0,0,0,0}
typealias Luminosity{T} SIQuantity{T,2,1,-3,0,0,0,0,0,0}
typealias Flux{T} SIQuantity{T,0,1,-3,0,0,0,0,0,0}
typealias Pressure{T} SIQuantity{T,-1,1,-2,0,0,0,0,0,0}
typealias Opacity{T} SIQuantity{T,2,-1,0,0,0,0,0,0,0}
typealias Acceleration{T} SIQuantity{T,1,0,-2,0,0,0,0,0,0}
typealias Density{T} SIQuantity{T,-3,1,0,0,0,0,0,0,0}
typealias Gravity Acceleration
typealias Dimensionless Real
typealias OpticalDepth Dimensionless
=#

# no explicit unit checking at the moment
typealias Mass Real
typealias Distance Real
typealias Radius Distance
typealias Temperature Real
typealias MassLuminosity Real
typealias Luminosity Real
typealias Flux Real
typealias Pressure Real
typealias Opacity Real
typealias Acceleration Real
typealias Density Real
typealias Gravity Acceleration
typealias Dimensionless Real
typealias OpticalDepth Dimensionless
