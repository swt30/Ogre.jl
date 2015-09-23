# The heat capacity of materials, as required for the adiabatic temperature gradient

using Dierckx


# Types for different heat capacity behaviours

"An isobaric heat capacity (cₚ)"
abstract HeatCapacity

"A constant heat capacity (no variation with T or P)"
immutable ConstantHeatCapacity{T<:Real} <: HeatCapacity
    value::T
end

"A heat capacity that's not constant"
abstract VaryingHeatCapacity <: HeatCapacity
abstract FunctionalHeatCapacity <: VaryingHeatCapacity

"A heat capacity which is a function of temperature"
immutable TFuncHeatCapacity{F} <: FunctionalHeatCapacity
    func::F
end

"A heat capacity which is a function of pressure and temperature"
immutable PTFuncHeatCapacity{F} <: FunctionalHeatCapacity
    func::F
end

"A heat capacity interpolated from a log-linear grid"
immutable GridHeatCapacity <: HeatCapacity
    P::Vector{Float64}
    T::Vector{Float64}
    spline::Spline2D

    function GridHeatCapacity(P, T, Cₚ)
        new(P, T, Spline2D(P, T, Cₚ, kx=1, ky=1))
    end
end


# Making heat capacities from files

"Generate a heat capacity from a file"
function GridHeatCapacity(filename)
    data = readdlm(filename, Float64)
    T = vec(data[2:end, 1])
    P = vec(data[1, 2:end])
    Cₚ = Matrix{Float64}(data[2:end, 2:end])' .* 1000  # note the transpose

    GridHeatCapacity(P, T, Cₚ)
end


# Evaluating heat capacities

Base.call(cp::GridHeatCapacity, P::Number, T::Number) = evaluate(cp.spline, P, T)
Base.call(cp::HeatCapacity, pv::PhysicalValues) = cp(pressure(pv), temperature(pv))
Base.call(cp::ConstantHeatCapacity, P::Number) = cp.value
Base.call(cp::ConstantHeatCapacity, P::Number, T::Number) = cp.value
Base.call(cp::TFuncHeatCapacity, T::Number) = cp.func(T)
Base.call(cp::TFuncHeatCapacity, P::Number, T::Number) = cp.func(T)
Base.call(cp::PTFuncHeatCapacity, P::Number, T::Number) = cp.func(P, T)
