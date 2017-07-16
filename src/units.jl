# Definitions of physical unit types for dispatch. At the moment this is a
# "dummy" set of definitions that just sets SI units to 1, but this may be
# swapped out for actual units if I have time.

using BasicUnits

const Mass = Real
const Distance = Real
const Radius = Distance
const Temperature = Real
const MassLuminosity = Real
const Luminosity = Real
const Flux = Real
const Pressure = Real
const Opacity = Real
const Acceleration = Real
const Density = Real
const Gravity = Acceleration
const Dimensionless = Real
const OpticalDepth = Dimensionless

#= definitions unused for now
using SIUnits
using SIUnits.ShortUnits
import SIUnits.unit, SIUnits.SIQuantity

# const Unit{T} = SIQuantity{T,...}
const Mass{T} = SIQuantity{T,0,1,0,0,0,0,0,0,0}
const Distance{T} = SIQuantity{T,1,0,0,0,0,0,0,0,0}
const Radius = Distance
const Temperature{T} = SIQuantity{T,0,0,0,0,1,0,0,0,0}
const MassLuminosity{T} = SIQuantity{T,2,0,-3,0,0,0,0,0,0}
const Luminosity{T} = SIQuantity{T,2,1,-3,0,0,0,0,0,0}
const Flux{T} = SIQuantity{T,0,1,-3,0,0,0,0,0,0}
const Pressure{T} = SIQuantity{T,-1,1,-2,0,0,0,0,0,0}
const Opacity{T} = SIQuantity{T,2,-1,0,0,0,0,0,0,0}
const Acceleration{T} = SIQuantity{T,1,0,-2,0,0,0,0,0,0}
const Density{T} = SIQuantity{T,-3,1,0,0,0,0,0,0,0}
const Gravity = Acceleration
const Dimensionless = Real
const OpticalDepth = Dimensionless
=#
