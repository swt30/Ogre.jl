using Base.Test
import Ogre
import WaterData


module test_structure_resources
import Ogre
pressure = Ogre.pressure
temperature = Ogre.temperature

type PressureEOS <: Ogre.EOS end
Ogre.@addEOSCall PressureEOS
(::PressureEOS)(P) = 4100. + 0.00161*(P^0.541)
Ogre.istempdependent(::PressureEOS) = false

type PTEOS <: Ogre.EOS end
Ogre.@addEOSCall PTEOS
(::PTEOS)(P, T) = 4100. + 0.00161*(P^0.541)
# the effect of T on ρ here is arbitrary - we just want increasing
# T to mean decreasing ρ
Ogre.istempdependent(::PTEOS) = true
end


@testset "Structure equations" begin
    res = test_structure_resources

    # realistic values taken from PREM at r=1000 km from the centre
    mass = 5.4e22      # kg enclosed
    radius = 1e6       # m
    pressure = 340e9   # Pa

    @testset "with no temperature dependence" begin
        realistic_values = Ogre.ValueSet(mass, radius, pressure)

        @testset "using an analytic EOS" begin
            eos = res.PressureEOS()

            masscontinuity = Ogre.MassContinuity(eos)
            pressurebalance = Ogre.PressureBalance()

            @testset "Zero gradients at the r=0 limit" begin
                zero_values = zero(Ogre.MassRadiusPressure)
                @test pressurebalance(zero_values) == 0
                @test masscontinuity(zero_values) == 0
            end

            @testset "Correct signs for the structure equations" begin
                # change in pressure is negative outwards
                @test pressurebalance(realistic_values) < 0
                # change in mass is positive outwards
                @test masscontinuity(realistic_values) > 0
            end

            @testset "Attempting to go over the singularity at r=0 returns zero" begin
                negative_radius = Ogre.ValueSet(5.4e22, -1e6, 340e9)
                negative_pressure = Ogre.ValueSet(5.4e22, 1e6, -340e9)
                negative_mass = Ogre.ValueSet(-5.4e22, 1e6, 340e9)
                map([negative_radius, negative_pressure, negative_mass]) do vs
                    @test pressurebalance(vs) == 0
                    @test masscontinuity(vs) == 0
                end
            end
        end
    end

    @testset "with temperature dependence" begin
        # realistic values taken from PREM at r=1000 km from the centre
        mass = 5.4e22      # kg enclosed
        radius = 1e6       # m
        pressure = 340e9   # Pa
        temperature = 5000 # K, this value is only rough
        realistic_values = Ogre.ValueSet(mass, radius, pressure, temperature)

        @testset "using an analytic EOS" begin
            eos = res.PTEOS()

            # arbitrary heat capacity
            Cₚ = WaterData.ConstantHeatCapacity(100)

            masscontinuity = Ogre.MassContinuity(eos)
            pressurebalance = Ogre.PressureBalance()
            temperaturegradient = Ogre.TemperatureGradient(eos, Cₚ)

            @testset "Zero gradients at the r=0 limit" begin
                zero_values = zero(Ogre.PhysicalValues)
                @test pressurebalance(zero_values) == 0
                @test masscontinuity(zero_values) == 0
                @test temperaturegradient(zero_values) == 0
            end

            @testset "Correct signs for the structure equations" begin
                # change in pressure is negative outwards
                @test pressurebalance(realistic_values) < 0
                # change in mass is positive outwards
                @test masscontinuity(realistic_values) > 0
                # change in temperature is negative outwards
                @test temperaturegradient(realistic_values) < 0
            end

            @testset "Attempting to go over the singularity at r=0 returns zero" begin
                negative_radius = Ogre.ValueSet(5.4e22, -1e6, 340e9, 5000.)
                negative_pressure = Ogre.ValueSet(5.4e22, 1e6, -340e9, 5000.)
                negative_mass = Ogre.ValueSet(-5.4e22, 1e6, 340e9, 5000.)
                negative_temperature = Ogre.ValueSet(5.4e22, -1e6, 340e9, -5000.)
                map([negative_radius, negative_pressure,
                     negative_mass, negative_temperature]) do vs
                         @test pressurebalance(vs) == 0
                         @test masscontinuity(vs) == 0
                         @test temperaturegradient(vs) == 0
                     end
            end
        end
    end
end


@testset "Planetary structure types" begin
    res = test_structure_resources
    @testset "Boundary values" begin
        bvs = Ogre.BoundaryValues(1,2,3)
        bvsT = Ogre.BoundaryValues(1,2,3,4)
        @test Ogre.radius(bvs) == 2
        @test Ogre.temperature(bvsT) == 4
    end

    @testset "Centre and surface values" begin
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

        @test Ogre.mass(Ogre.centre(s1)) == mc
        @test Ogre.mass(Ogre.centre(s2)) == mc
        @test Ogre.mass(Ogre.surface(s1)) == ms
        @test Ogre.mass(Ogre.surface(s2)) == ms
        @test Ogre.radius(Ogre.centre(s1)) == rc
        @test Ogre.radius(Ogre.centre(s2)) == rc
        @test Ogre.radius(Ogre.surface(s1)) == rs
        @test Ogre.radius(Ogre.surface(s2)) == rs
        @test Ogre.pressure(Ogre.centre(s1)) == Pc
        @test Ogre.pressure(Ogre.centre(s2)) == Pc
        @test Ogre.pressure(Ogre.surface(s1)) == Ps
        @test Ogre.pressure(Ogre.surface(s2)) == Ps
        @test Ogre.temperature(Ogre.centre(s2)) == Tc
        @test Ogre.temperature(Ogre.surface(s2)) == Ts
    end

    @testset "Planet system and solution setup" begin
        M = 5.972e24
        R = 6.3781e6
        Psurf = 1e5
        solution_grid = linspace(M, 0, 5)
        radius_bracket = [0, 10] * M

        pressurebalance = Ogre.PressureBalance()

        @testset "No temperature dependence" begin
            bvs = Ogre.BoundaryValues(M, R, Psurf)
            eos = res.PressureEOS()
            masscontinuity = Ogre.MassContinuity(eos)
            structure = Ogre.EquationSet([masscontinuity,
                                          pressurebalance])
            system = Ogre.PlanetSystem(M, eos, bvs, solution_grid,
                                       radius_bracket)
            struct = Ogre.blank_structure(system)

            @test length(Ogre.mass(struct)) == 5
            @test length(Ogre.radius(struct)) == 5
            @test length(Ogre.pressure(struct)) == 5
            @test_throws MethodError Ogre.temperature(struct)
        end

        @testset "Temperature dependence" begin
            Tsurf = 300
            bvs = Ogre.BoundaryValues(M, R, Psurf, Tsurf)
            eos = res.PTEOS()
            Cₚ = WaterData.ConstantHeatCapacity(1000)
            masscontinuity = Ogre.MassContinuity(eos)
            temperaturegradient = Ogre.TemperatureGradient(eos, Cₚ)
            structure = Ogre.EquationSet([masscontinuity,
                                          pressurebalance,
                                          temperaturegradient])
            system = Ogre.PlanetSystem(M, eos, Cₚ, bvs, solution_grid,
                                       radius_bracket)
            struct = Ogre.blank_structure(system)

            @test length(Ogre.mass(struct)) == 5
            @test length(Ogre.radius(struct)) == 5
            @test length(Ogre.pressure(struct)) == 5
            @test length(Ogre.temperature(struct)) == 5
        end
    end
end
