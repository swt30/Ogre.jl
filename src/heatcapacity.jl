# HEATCAPACITY.JL
# Types related to the heat capacity of materials, as required for the
# adiabatic temperature gradient

type HeatCapacity <: Equation
    equation::Function
end

Base.call(cp::HeatCapacity, x::Real) = cp.equation(x)
Base.call(cp::HeatCapacity, pv::PhysicalValues) = cp(pv.T)
