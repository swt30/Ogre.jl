using Ogre: Common, Structure, Eos
using FactCheck

facts("Structure equations") do
    # realistic values taken from PREM at r=1000 km
    realistic_values = ValueSet(5.4e22, 1e6, 340e9)
    eos_func(P::Real) = 4100. + 0.00161*(P^0.541)
    eos = Eos.SimpleEOS(eos_func, "Analytic MgSiO3")

    context("Zero gradients at the r=0 limit") do
        zero_values = zero(ValueSet)
        @fact Structure.pressure_balance(zero_values) => 0
        @fact Structure.mass_continuity(zero_values, eos) => 0
    end

    context("Correct signs for the structure equations") do
        # change in pressure is negative outwards
        @fact Structure.pressure_balance(realistic_values) => less_than(0)
        # change in mass is positive outwards
        @fact Structure.mass_continuity(realistic_values, eos) => greater_than(0)
    end

    context("Error if we attempt to go over the singularity at r=0") do
        negative_radius = ValueSet(5.4e22, -1e6, 340e9)
        negative_pressure = ValueSet(5.4e22, 1e6, -340e9)
        negative_mass = ValueSet(-5.4e22, 1e6, 340e9)
        map([negative_radius, negative_pressure, negative_mass]) do vs
            @fact_throws Structure.pressure_balance(vs)
            @fact_throws Structure.mass_continuity(vs, eos)
        end
    end
end
