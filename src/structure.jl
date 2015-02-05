module Structure
using Ogre: Common, Constants, Eos
# Exported types
export StructureEquation
# Exported functions
export mass_continuity, pressure_balance

immutable StructureEquation <: Equation
    equation::Function
end

isnegative(val::Real) = val < 0 ? true : false

function mass_continuity(vs::ValueSet, eos::EOS)
    m, r, P = vs.m, vs.r, vs.P

    if m <= 0 || r <= 0 || P <= 0
        return 0.0
    end

    rho = callfunc(eos, vs)
    dr_dm = 1 / (4 * pi * r^2 * rho)
    dr_dm::Float64
end

function pressure_balance(vs::ValueSet)
    m, r, P = vs.m, vs.r, vs.P

    if m <= 0 || r <= 0 || P <= 0
        return 0.0
    end

    dP_dm = -(G * m) / (4 * pi * r^4)
    dP_dm::Float64
end

end
