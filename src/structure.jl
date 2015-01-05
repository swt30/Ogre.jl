module Structure
using Ogre: Common, Constants, Eos
# Exported types
export StructureEquation
# Exported functions
export mass_continuity, pressure_balance

immutable StructureEquation <: Equation
    equation::Function
end

function mass_continuity(vs::ValueSet, eos::EOS)
    m, r, P = vs.m, vs.r, vs.P
    rho = callfunc(eos, vs)

    dr_dm = 0.

    if r > 0
        dr_dm = 1 / (4 * pi * r^2 * rho)
    end

    dr_dm::Float64
end

function pressure_balance(vs::ValueSet)
    m, r, P = vs.m, vs.r, vs.P

    dP_dm = 0.

    if r > 0
        dP_dm = -(G * m) / (4 * pi * r^4)
    end

    dP_dm::Float64
end

end
