@doc "Type for a planetary structure equation that is not an EOS" ->
immutable StructureEquation <: Equation
    equation::Function
end

@doc "The mass continuity equation: dr/dm = 1/4πr²ρ" ->
function mass_continuity{T<:Real}(vs::ValueSet{T}, eos::EOS)
    if vs.m <= 0 || vs.r <= 0 || vs.P <= 0
        return 0.0
    end

    rho = eos(vs)
    dr_dm::Float64 = 1 / (4 * pi * vs.r^2 * rho)
end

@doc "The pressure balance equation: dP/dm = -Gm/4πr⁴" ->
function pressure_balance{T<:Real}(vs::ValueSet{T})
    if vs.m <= 0 || vs.r <= 0 || vs.P <= 0
        return 0.0
    end

    dP_dm::Float64 = -(G * vs.m) / (4 * pi * vs.r^4)
end
