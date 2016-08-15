using Base.Test
import Ogre
import WaterData

module test_integrator_resources
import Ogre
import WaterData

# Simple ODEs
# dy/dt = 6  -->  y = 6t (y0=0)
ode_1(t, y::Real) = 6
ode_1(t, y::Vector) = [6]
ode_1_sol(t) = 6.*t
# dy/dt = 2t -->  t = t^2 + y0 (y0=0)
ode_2(t, y) = 2.*t
ode_2_sol(t) = t.^2
# dy/dt = y  -->  y = y0 e^t (y0=1)
ode_3(t, y) = y
ode_3_sol(t) = exp(t)
# dy1/dt = -y2, dy2/dt = y1  --> oscillating solution in y and w
ode_4(t, y) = [-y[2], y[1]]
ode_4_sol(t) = hcat(cos(t)-2*sin(t), 2*cos(t) + sin(t))

simple_odes = [ode_1, ode_2, ode_3, ode_4]
simple_ode_solutions = [ode_1_sol, ode_2_sol, ode_3_sol, ode_4_sol]

# Simple no-temperature EOSes
fe = let
    vinet = Ogre.Vinet(8.30e3, 156.2e9, 6.08, [1e3, 1e6]...)
    tfd = Ogre.TFD(26, 55.845)
    Ogre.PressurePiecewiseEOS([vinet, tfd], [0, 2.09e13, Inf])
end
mgsio3 = let
    bme = Ogre.BME3(4.10e3, 247e9, 3.97, [1e3, 1e6]...)
    tfd = Ogre.TFD([12, 14, 8], [24.305, 28.0855, 15.9994], [1., 1., 3.])
    Ogre.PressurePiecewiseEOS([bme, tfd], [0, 1.35e13, Inf])
end
h2o = let
    bme = Ogre.BME3(1460., 23.7e9, 4.15, [1e3, 1e6]...)
    tfd = Ogre.TFD([1, 8], [1.00794, 15.9994], [2., 1.])
    Ogre.PressurePiecewiseEOS([bme, tfd], [0, 7686e9, Inf])
end

# Simple with-temperature EOS
type JustAPressureEOS <: Ogre.EOS end
Ogre.@addEOSCall JustAPressureEOS
PV = Ogre.PhysicalValues
MRP = Ogre.MassRadiusPressure
Base.call(eos::JustAPressureEOS, pv::PV) = eos(Ogre.pressure(pv))
Base.call(eos::JustAPressureEOS, pv::MRP) = eos(Ogre.pressure(pv))
Base.call(::JustAPressureEOS, P) = 1000 + P
Ogre.istempdependent(::JustAPressureEOS) = false

type NotReallyATemperatureEOS <: Ogre.EOS end
@Ogre.addEOSCall NotReallyATemperatureEOS
Base.call(::NotReallyATemperatureEOS, P, T) = 1000 + P
Ogre.istempdependent(::NotReallyATemperatureEOS) = true

type ReallyATemperatureEOS <: Ogre.EOS end
@Ogre.addEOSCall ReallyATemperatureEOS
Base.call(::ReallyATemperatureEOS, P, T) = 1000 + 1/T
Ogre.istempdependent(::ReallyATemperatureEOS) = true

# Planet structures
M_earth = Ogre.M_earth
R_earth = Ogre.R_earth
atmospheric_pressure = 1e5
surface_temperature = 300
solution_grid = linspace(M_earth, 0, 100)
r_bracket = [0, 2]*M_earth
dual_layer_eos = let
    eoses = [mgsio3, fe]
    transitions = [0, 2/3, 1] * M_earth
    Ogre.MassPiecewiseEOS(eoses, transitions)
end
tri_layer_eos = let
    eoses = [h2o, mgsio3, fe]
    transitions = [0, 1/6, 2/3, 1] * M_earth
    Ogre.MassPiecewiseEOS(eoses, transitions)
end
tri_layer_system = let
    trial_surface = Ogre.BoundaryValues(M_earth, R_earth, atmospheric_pressure)
    system = Ogre.PlanetSystem(M_earth, tri_layer_eos, trial_surface,
                               solution_grid)
    Ogre.converge!(system)
    system
end
tri_layer_radius = Ogre.currentradiusguess(tri_layer_system)

# Full temperature structure
withtemp = let
    heatcap = WaterData.ConstantHeatCapacity(4200)
    eos = ReallyATemperatureEOS()
    surface = Ogre.BoundaryValues(M_earth, R_earth, atmospheric_pressure,
                                  surface_temperature)
    system = Ogre.PlanetSystem(M_earth, eos, heatcap, surface, solution_grid,
                               r_bracket)
    Ogre.converge!(system)
    system
end
end


@testset "Integrator tests" begin
    mass = Ogre.mass
    radius = Ogre.radius
    pressure = Ogre.pressure
    temperature = Ogre.temperature
    centre = Ogre.centre
    surface = Ogre.surface
    res = test_integrator_resources

    @testset "Standard RK4 integrator test cases" begin
        t = collect(Float64, 0:.1:10)

        # Four different integrator tests
        odes = res.simple_odes
        soln_funcs = res.simple_ode_solutions
        solns = map(f -> f(t), soln_funcs)

        @testset "The integrator initialises and steps properly" begin
            integrator = Ogre.RK4
            s1 = integrator(odes[1], 0., t)
            s2 = integrator(odes[1], [0.], t)
            s3 = integrator(odes[4], [1., 2.], t)
            @test Ogre.y0(s1) == 0.0
            @test start(s1) == 1
            @test next(s1, 1)[2] == 2
            @test Ogre.y0(s2) == [0.0]
            @test start(s2) == 1
            @test next(s2, 1)[2] == 2
            @test Ogre.y0(s3) == [1., 2.]
            @test start(s3) == 1
            @test next(s3, 1)[2] == 2
        end

        @testset "The dense solutions are as expected" begin
            f = Ogre.ode4_dense
            @test isapprox(f(odes[1], 0., t), solns[1])
            @test isapprox(f(odes[2], 0., t), solns[2])
            @test isapprox(f(odes[3], 1., t), solns[3], rtol=5e-3)
            @test isapprox(f(odes[4], [1., 2.], t), solns[4], rtol=5e-3)
        end
    end

    @testset "Planetary integrator test cases" begin
        @testset "Solving for a radius (high-level funcs)" begin
            @testset "We can match Madhu and Sara's results to within 1%" begin
                test_masses = [0.5, 1.0, 5.0]
                fe = res.fe
                mgsio3 = res.mgsio3
                target_fe_radii = [0.628896, 0.76727, 1.17721]
                target_mgsio3_radii = [0.842621, 1.04239, 1.65182]

                fe_radii = Ogre.R(test_masses, fe, in_earth_units=true)
                mgsio3_radii = Ogre.R(test_masses, mgsio3, in_earth_units=true)

                @test isapprox(fe_radii, target_fe_radii, rtol=0.01)
                @test isapprox(mgsio3_radii, target_mgsio3_radii, rtol=0.01)
            end

            @testset "Dual- and tri-layer solutions work" begin
                # accept a pretty large error since we're really just looking
                # to make sure that it produces the right type of solution
                @test isapprox(Ogre.R(res.M_earth, res.dual_layer_eos),
                               res.R_earth, rtol=0.05)
                @test isapprox(Ogre.R(res.M_earth, res.tri_layer_eos),
                               res.R_earth, rtol=0.05)
            end
        end

        @testset "Solving for internal structure (tri-layer)" begin
            system = res.tri_layer_system
            soln = Ogre.solve(system)

            @testset "Solution outputs are the correct type" begin
                @test isa(soln, Ogre.PlanetStructure{Ogre.NoTemp})
            end

            @testset "Solution boundaries match the boundary conditions" begin
                @test mass(centre(soln)) == 0
                @test mass(surface(soln)) == res.M_earth
                @test pressure(centre(soln)) > pressure(surface(soln))
                @test pressure(surface(soln)) == res.atmospheric_pressure
                @test radius(centre(soln)) < 100
                @test radius(centre(soln)) > 0
                @test radius(surface(soln)) == res.tri_layer_radius
            end
        end

        @testset "Solving a single-layer solution two different ways" begin
            eos_notemp = res.JustAPressureEOS()
            eos_withtemp = res.NotReallyATemperatureEOS()
            heatcap = WaterData.ConstantHeatCapacity(1000.)
            local radius = Ogre.R(res.M_earth, eos_notemp)
            tempdep_radius = Ogre.R(res.M_earth, eos_withtemp, heatcap)
            @test tempdep_radius == radius
        end

        @testset "Temperature dependence" begin
            system = res.withtemp
            soln = Ogre.solve(system)

            @test mass(centre(soln)) == 0
            @test mass(surface(soln)) == res.M_earth
            @test pressure(centre(soln)) > pressure(surface(soln))
            @test pressure(surface(soln)) == res.atmospheric_pressure
            @test temperature(centre(soln)) > temperature(surface(soln))
            @test temperature(surface(soln)) == 300
            @test radius(centre(soln)) < 100
            @test isa(soln, Ogre.PlanetStructure{Ogre.WithTemp})
        end
    end
end
