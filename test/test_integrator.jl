include("header.jl")

facts("Integrator tests") do

    mass = Ogre.mass
    radius = Ogre.radius
    pressure = Ogre.pressure
    temperature = Ogre.temperature
    centre = Ogre.centre
    surface = Ogre.surface

    context("Standard RK4 integrator test cases") do
        t = collect(0:0.1:10)

        # Four different integrator tests
        odes = res.simple_odes
        soln_funcs = res.simple_ode_solutions
        solns = map(f -> f(t), soln_funcs)

        context("The integrator initialises properly") do
            integrator = Ogre.GenericRK4
            s0 = integrator(odes[1], 0., t)
            s1 = integrator(odes[1], [0.], t)
            s2 = integrator(odes[4], [1., 2.], t)
            @fact start(s0) => (1, 0.0)
            @fact start(s1) => (1, [0.0])
            @fact start(s2) => (1, [1., 2.])
        end

        context("The dense solutions are as expected") do
            solver = Ogre.ode4_dense
            @fact solver(odes[1], 0., t) => roughly(solns[1])
            @fact solver(odes[2], 0., t) => roughly(solns[2])
            @fact solver(odes[3], 1., t) => roughly(solns[3])
            @fact solver(odes[4], [1., 2.], t) => roughly(solns[4], rtol = 1e-3)
        end
    end

    context("Planetary integrator test cases") do
        context("Solving for a radius (high-level funcs)") do
            context("We can match Madhu and Sara's results to within 1%") do
                test_masses = [0.5, 1.0, 5.0]
                fe = Ogre.fe_seager
                mgsio3 = Ogre.mgsio3_seager
                target_fe_radii = [0.628896, 0.76727, 1.17721]
                target_mgsio3_radii = [0.842621, 1.04239, 1.65182]

                fe_radii = Ogre.R(test_masses, fe, in_earth_units=true)
                mgsio3_radii = Ogre.R(test_masses, mgsio3, in_earth_units=true)

                @fact fe_radii => roughly(target_fe_radii, rtol=0.01)
                @fact mgsio3_radii => roughly(target_mgsio3_radii, rtol=0.01)
            end

            context("Dual- and tri-layer solutions work") do
                # accept a pretty large error since we're really just looking
                # to make sure that it produces the right type of solution
                @fact (Ogre.R(res.planet.M_earth, res.planet.dual_layer.eos)
                       => roughly(res.planet.R_earth, rtol=0.05))
                @fact (Ogre.R(res.planet.M_earth, res.planet.tri_layer.eos)
                       => roughly(res.planet.R_earth, rtol=0.05))
            end
        end

        context("Solving for internal structure (tri-layer)") do
            system = res.planet.tri_layer.system
            soln = Ogre.solve(system)

            context("Solution outputs are the correct type") do
                @fact isa(soln, Ogre.PlanetStructure{Ogre.NoTemp}) => true
            end

            context("Solution boundaries match the boundary conditions") do
                @fact mass(centre(soln)) => 0
                @fact mass(surface(soln)) => res.planet.M_earth
                @fact pressure(centre(soln)) => greater_than(pressure(surface(soln)))
                @fact pressure(surface(soln)) => res.planet.atmospheric_pressure
                @fact radius(centre(soln)) => less_than(100)
                @fact radius(centre(soln)) => greater_than(0)
                @fact radius(surface(soln)) => res.planet.tri_layer.actual_radius
            end
        end

        context("Solving a single-layer solution two different ways") do
            f(P) = Ogre.mgsio3_func(P)
            f(P, T) = f(P)
            eos_notemp = Ogre.SimpleEOS(Ogre.NoTemp, f, "No temperature dependence")
            eos_withtemp = Ogre.SimpleEOS(Ogre.WithTemp, f, "With temp dependence")
            actual_radius = Ogre.R(res.planet.M_earth, eos_notemp)
            tempdep_radius = Ogre.R(res.planet.M_earth, eos_withtemp, res.heatcap.exponential)
            @fact tempdep_radius => actual_radius
        end

        context("Helper function tests") do
            system = Ogre.DefaultPlanetSystem(res.planet.M_earth, res.planet.dual_layer.eos)

            context("Determining when we're close to the centre of the planet") do
                @fact Ogre.hit_the_centre(-50) => true
                @fact Ogre.hit_the_centre(50) => false
                @fact Ogre.hit_the_centre(1e6) => false
                @fact Ogre.not_far_enough(1e6) => true
                @fact Ogre.not_far_enough(50) => false
            end
        end

        context("Temperature dependence") do
            system = res.planet.withtemp.system
            soln = Ogre.solve(system)

            @fact mass(centre(soln)) => 0
            @fact mass(surface(soln)) => res.planet.M_earth
            @fact pressure(centre(soln)) => greater_than(pressure(surface(soln)))
            @fact pressure(surface(soln)) => res.planet.atmospheric_pressure
            @fact temperature(centre(soln)) => greater_than(temperature(surface(soln)))
            @fact temperature(surface(soln)) => 300
            @fact radius(centre(soln)) => less_than(100)
            @fact isa(soln, Ogre.PlanetStructure{Ogre.WithTemp}) => true
        end
    end
end


