# HEATCAPACITY.JL
# Types related to the heat capacity of materials, as required for the
# adiabatic temperature gradient

using Dierckx

@doc "An isobaric heat capacity (Cₚ)" ->
type HeatCapacity <: Equation
    equation::Function
end
HeatCapacity(::Type{WithTemp}, f::Function) = HeatCapacity(f)

import Base.call
call(cp::HeatCapacity, x::Real) = cp.equation(x)
call(cp::HeatCapacity, pv::PhysicalValues) = cp(temperature(pv))

function HeatCapacity(filename::String)
    data = readdlm(filename)
    T = Vector{Float64}(data[2:end, 1])
    P = Vector{Float64}(vec(data[1, 2:end]))
    Cₚs = Matrix{Float64}(data[2:end, 2:end]) .* 1000

    Cₚ = Cₚs[:, 6] # 1: low temp, 6: hi-temp
    heatcap_f = lininterp(T, Cₚ)

    HeatCapacity(heatcap_f)
end