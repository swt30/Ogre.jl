import Ogre
using FactCheck

facts("Structure equations") do
    # realistic values taken from PREM at r=1000 km from the centre
    mass = 5.4e22      # kg enclosed
    radius = 1e6       # m
    pressure = 340e9   # Pa

    context("with no temperature dependence") do
        realistic_values = Ogre.ValueSet(mass, radius, pressure)

        context("using an analytic EOS") do
            eos_func(P::Real) = 4100. + 0.00161*(P^0.541)
            eos = Ogre.PressureEOS(eos_func, "Analytic MgSiO3")

            context("Zero gradients at the r=0 limit") do
                zero_values = zero(Ogre.MassRadiusPressure)
                @fact Ogre.pressure_balance(zero_values) => 0
                @fact Ogre.mass_continuity(zero_values, eos) => 0
            end

            context("Correct signs for the structure equations") do
                # change in pressure is negative outwards
                @fact Ogre.pressure_balance(realistic_values) => less_than(0)
                # change in mass is positive outwards
                @fact Ogre.mass_continuity(realistic_values, eos) => greater_than(0)
            end

            context("Attempting to go over the singularity at r=0 returns zero") do
                negative_radius = Ogre.ValueSet(5.4e22, -1e6, 340e9)
                negative_pressure = Ogre.ValueSet(5.4e22, 1e6, -340e9)
                negative_mass = Ogre.ValueSet(-5.4e22, 1e6, 340e9)
                map([negative_radius, negative_pressure, negative_mass]) do vs
                    @fact Ogre.pressure_balance(vs) => 0
                    @fact Ogre.mass_continuity(vs, eos) => 0
                end
            end
        end
    end

    context("with temperature dependence") do
        # realistic values taken from PREM at r=1000 km from the centre
        mass = 5.4e22      # kg enclosed
        radius = 1e6       # m
        pressure = 340e9   # Pa
        temperature = 5000 # K, this value is only rough
        realistic_values = Ogre.ValueSet(mass, radius, pressure, temperature)

        context("using an analytic EOS") do
            eos_func(P::Real, T::Real) = 4100. + 0.00161*(P^0.541)*(T^-0.054)
            eos = Ogre.PressureTempEOS(eos_func, "Analytic MgSiO3 with temp")
            # the effect of T on ρ here is arbitrary - we just want increasing
            # T to mean decreasing ρ

            # arbitrary heat capacity
            heatcap(args...) = 100
            Cₚ = Ogre.HeatCapacity(heatcap)

            context("Zero gradients at the r=0 limit") do
                zero_values = zero(Ogre.PhysicalValues)
                @fact Ogre.pressure_balance(zero_values) => 0
                @fact Ogre.mass_continuity(zero_values, eos) => 0
                @fact Ogre.temperature_gradient(zero_values, eos, Cₚ) => 0
            end

            context("Correct signs for the structure equations") do
                # change in pressure is negative outwards
                @fact Ogre.pressure_balance(realistic_values) => less_than(0)
                # change in mass is positive outwards
                @fact Ogre.mass_continuity(realistic_values, eos) => greater_than(0)
                # change in temperature is negative outwards
                @fact Ogre.temperature_gradient(realistic_values, eos, Cₚ) => less_than(0)
            end

            context("Attempting to go over the singularity at r=0 returns zero") do
                negative_radius = Ogre.ValueSet(5.4e22, -1e6, 340e9, 5000.)
                negative_pressure = Ogre.ValueSet(5.4e22, 1e6, -340e9, 5000.)
                negative_mass = Ogre.ValueSet(-5.4e22, 1e6, 340e9, 5000.)
                negative_temperature = Ogre.ValueSet(5.4e22, -1e6, 340e9, -5000.)
                map([negative_radius, negative_pressure,
                     negative_mass, negative_temperature]) do vs
                    @fact Ogre.pressure_balance(vs) => 0
                    @fact Ogre.mass_continuity(vs, eos) => 0
                    @fact Ogre.temperature_gradient(vs, eos, Cₚ) => 0
                end
            end
        end
    end
end

facts("Planetary structure types") do
    @pending "test that we can initialise the types" => true
    # etc etc
end
