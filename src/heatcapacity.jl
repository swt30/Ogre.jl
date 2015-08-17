# HEATCAPACITY.JL
# Types related to the heat capacity of materials, as required for the
# adiabatic temperature gradient

using Dierckx

"An isobaric heat capacity (cₚ)"
abstract HeatCapacity

"A constant heat capacity (no variation with T or P)"
immutable ConstantHeatCapacity{T<:Real} <: HeatCapacity
    value::T
end

"A heat capacity that's not constant"
abstract VaryingHeatCapacity <: HeatCapacity
abstract FunctionalHeatCapacity <: VaryingHeatCapacity

"A heat capacity that varies with temperature"
immutable TFuncHeatCapacity{F} <: FunctionalHeatCapacity
    func::F
end

"A heat capacity that varies with pressure and temperature"
immutable PTFuncHeatCapacity{F} <: FunctionalHeatCapacity
    func::F
end

"A heat capacity interpolated from a log-linear grid"
immutable PTGridHeatCapacity <: HeatCapacity
    logP::Vector{Float64}
    T::Vector{Float64}
    spline::Spline2D

    function PTGridHeatCapacity(P, T, Cₚ)
        new(log10(P), T, Spline2D(log10(P), T, Cₚ, kx=1, ky=1))
    end
end
Base.call(cp::PTGridHeatCapacity, P, T) = evaluate(cp.spline, log10(P), T)

HeatCapacity(::Type{WithTemp}, func) = TFuncHeatCapacity(func)
HeatCapacity(::Type{WithTempPressure}, func) = PTFuncHeatCapacity(func)
HeatCapacity(n::Number) = ConstantHeatCapacity(n)
function HeatCapacity(filename::String)
    # generate heat capacity from file
    # TODO: document this constructor once permitted
    data = readdlm(filename, Float64)
    T = vec(data[2:end, 1])
    P = vec(data[1, 2:end])
    Cₚ = Matrix{Float64}(data[2:end, 2:end])' .* 1000  # note the transpose

    PTGridHeatCapacity(P, T, Cₚ)
end

import Base.call
call(cp::HeatCapacity, pv::PhysicalValues) = cp(pressure(pv), temperature(pv))
call(cp::ConstantHeatCapacity, T::Real) = cp.value
call(cp::ConstantHeatCapacity, P::Real, T::Real) = cp.value
call(cp::TFuncHeatCapacity, T::Real) = cp.func(T)
call(cp::TFuncHeatCapacity, P::Real, T::Real) = cp.func(T)
call(cp::PTFuncHeatCapacity, P::Real, T::Real) = cp.func(P, T)