include("header.jl")

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

            masscontinuity = Ogre.MassContinuityEq(eos)
            pressurebalance = Ogre.PressureBalanceEq()

            context("Zero gradients at the r=0 limit") do
                zero_values = zero(Ogre.MassRadiusPressure)
                @fact pressurebalance(zero_values) => 0
                @fact masscontinuity(zero_values) => 0
            end

            context("Correct signs for the structure equations") do
                # change in pressure is negative outwards
                @fact pressurebalance(realistic_values) => less_than(0)
                # change in mass is positive outwards
                @fact masscontinuity(realistic_values) => greater_than(0)
            end

            context("Attempting to go over the singularity at r=0 returns zero") do
                negative_radius = Ogre.ValueSet(5.4e22, -1e6, 340e9)
                negative_pressure = Ogre.ValueSet(5.4e22, 1e6, -340e9)
                negative_mass = Ogre.ValueSet(-5.4e22, 1e6, 340e9)
                map([negative_radius, negative_pressure, negative_mass]) do vs
                    @fact pressurebalance(vs) => 0
                    @fact masscontinuity(vs) => 0
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

            masscontinuity = Ogre.MassContinuityEq(eos)
            pressurebalance = Ogre.PressureBalanceEq()
            temperaturegradient = Ogre.TemperatureGradientEq(eos, Cₚ)

            context("Zero gradients at the r=0 limit") do
                zero_values = zero(Ogre.PhysicalValues)
                @fact pressurebalance(zero_values) => 0
                @fact masscontinuity(zero_values) => 0
                @fact temperaturegradient(zero_values) => 0
            end

            context("Correct signs for the structure equations") do
                # change in pressure is negative outwards
                @fact pressurebalance(realistic_values) => less_than(0)
                # change in mass is positive outwards
                @fact masscontinuity(realistic_values) => greater_than(0)
                # change in temperature is negative outwards
                @fact temperaturegradient(realistic_values) => less_than(0)
            end

            context("Attempting to go over the singularity at r=0 returns zero") do
                negative_radius = Ogre.ValueSet(5.4e22, -1e6, 340e9, 5000.)
                negative_pressure = Ogre.ValueSet(5.4e22, 1e6, -340e9, 5000.)
                negative_mass = Ogre.ValueSet(-5.4e22, 1e6, 340e9, 5000.)
                negative_temperature = Ogre.ValueSet(5.4e22, -1e6, 340e9, -5000.)
                map([negative_radius, negative_pressure,
                     negative_mass, negative_temperature]) do vs
                    @fact pressurebalance(vs) => 0
                    @fact masscontinuity(vs) => 0
                    @fact temperaturegradient(vs) => 0
                end
            end
        end
    end
end

facts("Planetary structure types") do
    context("Boundary values") do
        bvs = Ogre.BoundaryValues(1,2,3)
        bvsT = Ogre.BoundaryValues(1,2,3,4)
        @fact Ogre.radius(bvs) => 2
        @fact Ogre.temperature(bvsT) => 4
    end

    context("Centre and surface values") do
        mc, ms = 1, 2
        rc, rs = 3, 4
        Pc, Ps = 5, 6
        Tc, Ts = 7, 8
        npoints = 10

        m = linspace(ms, mc, npoints)
        r = linspace(rs, rc, npoints)
        P = linspace(Ps, Pc, npoints)
        T = linspace(Ts, Tc, npoints)

        s1 = Ogre.PlanetStructure(m, r, P)
        s2 = Ogre.PlanetStructure(m, r, P, T)

        @fact Ogre.mass(Ogre.centre(s1)) => mc
        @fact Ogre.mass(Ogre.centre(s2)) => mc
        @fact Ogre.mass(Ogre.surface(s1)) => ms
        @fact Ogre.mass(Ogre.surface(s2)) => ms
        @fact Ogre.radius(Ogre.centre(s1)) => rc
        @fact Ogre.radius(Ogre.centre(s2)) => rc
        @fact Ogre.radius(Ogre.surface(s1)) => rs
        @fact Ogre.radius(Ogre.surface(s2)) => rs
        @fact Ogre.pressure(Ogre.centre(s1)) => Pc
        @fact Ogre.pressure(Ogre.centre(s2)) => Pc
        @fact Ogre.pressure(Ogre.surface(s1)) => Ps
        @fact Ogre.pressure(Ogre.surface(s2)) => Ps
        @fact Ogre.temperature(Ogre.centre(s2)) => Tc
        @fact Ogre.temperature(Ogre.surface(s2)) => Ts
    end

    context("Planet system and solution setup") do
        M = 5.972e24
        R = 6.3781e6
        Psurf = 1e5
        solution_grid = linspace(M, 0, 5)
        radius_bracket = [0, 10] * M

        pressurebalance = Ogre.PressureBalanceEq()

        context("No temperature dependence") do
            bvs = Ogre.BoundaryValues(M, R, Psurf)
            eos_f(P) = 4100. + 0.00161*(P^0.541)
            eos = Ogre.SimpleEOS(Ogre.NoTemp, eos_f, "")
            masscontinuity = Ogre.MassContinuityEq(eos)
            structure = Ogre.EquationSet([masscontinuity,
                                          pressurebalance])
            system = Ogre.PlanetSystem(M, eos, bvs, solution_grid,
                                       radius_bracket)
            struct = Ogre.blank_structure(system)

            @fact length(Ogre.mass(struct)) => 5
            @fact length(Ogre.radius(struct)) => 5
            @fact length(Ogre.pressure(struct)) => 5
            @fact_throws Ogre.temperature(struct) => 5
        end

        context("Temperature dependence") do
            Tsurf = 300
            bvs = Ogre.BoundaryValues(M, R, Psurf, Tsurf)
            eos_f(P, T) = 4100. + 0.00161*(P^0.541)*(T^-0.054)
            eos = Ogre.SimpleEOS(Ogre.WithTemp, eos_f, "")
            Cₚ_f(T) = T
            Cₚ = Ogre.HeatCapacity(Cₚ_f)
            masscontinuity = Ogre.MassContinuityEq(eos)
            temperaturegradient = Ogre.TemperatureGradientEq(eos, Cₚ)
            structure = Ogre.EquationSet([masscontinuity,
                                         pressurebalance,
                                         temperaturegradient])
            system = Ogre.PlanetSystem(M, eos, Cₚ, bvs, solution_grid,
                                       radius_bracket)
            struct = Ogre.blank_structure(system)

            @fact length(Ogre.mass(struct)) => 5
            @fact length(Ogre.radius(struct)) => 5
            @fact length(Ogre.pressure(struct)) => 5
            @fact length(Ogre.temperature(struct)) => 5
        end
    end
end
