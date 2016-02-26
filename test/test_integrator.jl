using FactCheck
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
PV = Ogre.PhysicalValues
MRP = Ogre.MassRadiusPressure
Base.call(eos::JustAPressureEOS, pv::PV) = eos(Ogre.pressure(pv))
Base.call(eos::JustAPressureEOS, pv::MRP) = eos(Ogre.pressure(pv))
Base.call(::JustAPressureEOS, P) = 1000 + P
Ogre.istempdependent(::JustAPressureEOS) = false

type NotReallyATemperatureEOS <: Ogre.EOS end
Base.call(::NotReallyATemperatureEOS, P, T) = 1000 + P
Ogre.istempdependent(::NotReallyATemperatureEOS) = true

type ReallyATemperatureEOS <: Ogre.EOS end
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


facts("Integrator tests") do
    mass = Ogre.mass
    radius = Ogre.radius
    pressure = Ogre.pressure
    temperature = Ogre.temperature
    centre = Ogre.centre
    surface = Ogre.surface
    res = test_integrator_resources

    context("Standard RK4 integrator test cases") do
        t = collect(Float64, 0:.1:10)

        # Four different integrator tests
        odes = res.simple_odes
        soln_funcs = res.simple_ode_solutions
        solns = map(f -> f(t), soln_funcs)

        context("The integrator initialises and steps properly") do
            integrator = Ogre.RK4
            s1 = integrator(odes[1], 0., t)
            s2 = integrator(odes[1], [0.], t)
            s3 = integrator(odes[4], [1., 2.], t)
            @fact Ogre.y0(s1) --> 0.0
            @fact start(s1) --> 1
            @fact next(s1, 1)[2] --> 2
            @fact Ogre.y0(s2) --> [0.0]
            @fact start(s2) --> 1
            @fact next(s2, 1)[2] --> 2
            @fact Ogre.y0(s3) --> [1., 2.]
            @fact start(s3) --> 1
            @fact next(s3, 1)[2] --> 2
        end

        context("The dense solutions are as expected") do
            f = Ogre.ode4_dense
            @fact f(odes[1], 0., t) --> roughly(solns[1])
            @fact f(odes[2], 0., t) --> roughly(solns[2])
            @fact f(odes[3], 1., t) --> roughly(solns[3], rtol=5e-3)
            @fact f(odes[4], [1., 2.], t) --> roughly(solns[4], rtol=5e-3)
        end
    end

    context("Planetary integrator test cases") do
        context("Solving for a radius (high-level funcs)") do
            context("We can match Madhu and Sara's results to within 1%") do
                test_masses = [0.5, 1.0, 5.0]
                fe = res.fe
                mgsio3 = res.mgsio3
                target_fe_radii = [0.628896, 0.76727, 1.17721]
                target_mgsio3_radii = [0.842621, 1.04239, 1.65182]

                fe_radii = Ogre.R(test_masses, fe, in_earth_units=true)
                mgsio3_radii = Ogre.R(test_masses, mgsio3, in_earth_units=true)

                @fact fe_radii --> roughly(target_fe_radii, rtol=0.01)
                @fact mgsio3_radii --> roughly(target_mgsio3_radii, rtol=0.01)
            end

            context("Dual- and tri-layer solutions work") do
                # accept a pretty large error since we're really just looking
                # to make sure that it produces the right type of solution
                @fact (Ogre.R(res.M_earth, res.dual_layer_eos)
                       --> roughly(res.R_earth, rtol=0.05))
                @fact (Ogre.R(res.M_earth, res.tri_layer_eos)
                       --> roughly(res.R_earth, rtol=0.05))
            end
        end

        context("Solving for internal structure (tri-layer)") do
            system = res.tri_layer_system
            soln = Ogre.solve(system)

            context("Solution outputs are the correct type") do
                @fact isa(soln, Ogre.PlanetStructure{Ogre.NoTemp}) --> true
            end

            context("Solution boundaries match the boundary conditions") do
                @fact mass(centre(soln)) --> 0
                @fact mass(surface(soln)) --> res.M_earth
                @fact pressure(centre(soln)) --> greater_than(pressure(surface(soln)))
                @fact pressure(surface(soln)) --> res.atmospheric_pressure
                @fact radius(centre(soln)) --> less_than(100)
                @fact radius(centre(soln)) --> greater_than(0)
                @fact radius(surface(soln)) --> res.tri_layer_radius
            end
        end

        context("Solving a single-layer solution two different ways") do
            eos_notemp = res.JustAPressureEOS()
            eos_withtemp = res.NotReallyATemperatureEOS()
            heatcap = WaterData.ConstantHeatCapacity(1000.)
            local radius = Ogre.R(res.M_earth, eos_notemp)
            tempdep_radius = Ogre.R(res.M_earth, eos_withtemp, heatcap)
            @fact tempdep_radius --> radius
        end

        context("Temperature dependence") do
            system = res.withtemp
            soln = Ogre.solve(system)

            @fact mass(centre(soln)) --> 0
            @fact mass(surface(soln)) --> res.M_earth
            @fact pressure(centre(soln)) --> greater_than(pressure(surface(soln)))
            @fact pressure(surface(soln)) --> res.atmospheric_pressure
            @fact temperature(centre(soln)) --> greater_than(temperature(surface(soln)))
            @fact temperature(surface(soln)) --> 300
            @fact radius(centre(soln)) --> less_than(100)
            @fact isa(soln, Ogre.PlanetStructure{Ogre.WithTemp}) --> true
        end
    end
end
