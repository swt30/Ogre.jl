# HEATCAPACITY.JL
# Types related to the heat capacity of materials, as required for the
# adiabatic temperature gradient

using Dierckx

@doc "An isobaric heat capacity (Cₚ)" ->
abstract HeatCapacity

type ConstantHeatCapacity{T<:Real} <: HeatCapacity
    value::T
end

abstract VaryingHeatCapacity <: HeatCapacity

type THeatCapacity <: VaryingHeatCapacity
    func
end

type TPHeatCapacity <: VaryingHeatCapacity
    func
end

HeatCapacity(::Type{WithTemp}, f::Function) = THeatCapacity(f)
HeatCapacity(::Type{WithTempPressure}, f::Function) = TPHeatCapacity(f)
HeatCapacity(n::Number) = ConstantHeatCapacity(n)

import Base.call
call(cp::HeatCapacity, pv::PhysicalValues) = cp(temperature(pv), pressure(pv))
call(cp::ConstantHeatCapacity, T::Real) = cp.value
call(cp::ConstantHeatCapacity, T::Real, P::Real) = cp.value
call(cp::THeatCapacity, T::Real) = cp.func(T)
call(cp::THeatCapacity, T::Real, P::Real) = cp.func(T)
call(cp::TPHeatCapacity, T::Real, P::Real) = cp.func(T, P)

function HeatCapacity(filename::String)
    data = readdlm(filename)
    T = Vector{Float64}(data[2:end, 1])
    P = Vector{Float64}(vec(data[1, 2:end]))
    Cₚs = Matrix{Float64}(data[2:end, 2:end]) .* 1000

    heatcap_f = semiloginterpy(T, P, Cₚs)

    HeatCapacity(WithTempPressure, heatcap_f)
end
