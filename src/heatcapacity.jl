# HEATCAPACITY.JL
# Types related to the heat capacity of materials, as required for the
# adiabatic temperature gradient

type HeatCapacity <: Equation
    equation::Function
end
HeatCapacity(::Type{WithTemp}, f::Function) = HeatCapacity(f)

import Base.call
call(cp::HeatCapacity, x::Real) = cp.equation(x)
call(cp::HeatCapacity, pv::PhysicalValues) = cp(temperature(pv))
