# HEATCAPACITY.JL
# Types related to the heat capacity of materials, as required for the
# adiabatic temperature gradient

using Dierckx

"An isobaric heat capacity (cₚ)"
abstract HeatCapacity

"A constant heat capacity (no variation with T or P)"
type ConstantHeatCapacity{T<:Real} <: HeatCapacity
    value::T
end

"A heat capacity that's not constant"
abstract VaryingHeatCapacity <: HeatCapacity

"A heat capacity that varies with temperature"
type THeatCapacity <: VaryingHeatCapacity
    func
end

"A heat capacity that varies with pressure and temperature"
type PTHeatCapacity <: VaryingHeatCapacity
    func
end

HeatCapacity(::Type{WithTemp}, f::Function) = THeatCapacity(f)
HeatCapacity(::Type{WithTempPressure}, f::Function) = PTHeatCapacity(f)
HeatCapacity(n::Number) = ConstantHeatCapacity(n)
function HeatCapacity(filename::String)
    # generate heat capacity from file
    # TODO: document this constructor once permitted
    data = readdlm(filename)
    T = Vector{Float64}(data[2:end, 1])
    P = Vector{Float64}(vec(data[1, 2:end]))
    Cₚs = Matrix{Float64}(data[2:end, 2:end]) .* 1000

    heatcap_f = semiloginterpx(P, T, Cₚs')

    HeatCapacity(WithTempPressure, heatcap_f)
end

import Base.call
call(cp::HeatCapacity, pv::PhysicalValues) = cp(pressure(pv), temperature(pv))
call(cp::ConstantHeatCapacity, T::Real) = cp.value
call(cp::ConstantHeatCapacity, P::Real, T::Real) = cp.value
call(cp::THeatCapacity, T::Real) = cp.func(T)
call(cp::THeatCapacity, P::Real, T::Real) = cp.func(T)
call(cp::PTHeatCapacity, P::Real, T::Real) = cp.func(P, T)