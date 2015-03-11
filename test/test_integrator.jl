include("header.jl")

facts("Integrator tests") do

    context("Standard RK4 integrator test cases") do
        t = collect(0:0.1:10)

        # Four different integrator tests

        # dy/dt = 6  -->  y = 6t (y0=0)
        f1(t, y) = 6
        sol1 = 6.*t
        # dy/dt = 2t -->  t = t^2 (y0=0)
        f2(t, y) = 2t
        sol2 = t.^2
        # dy/dt = y  -->  y = y0 e^t (y0=1)
        f3(t, y) = y
        sol3 = exp(t)
        # dy1/dt = -y2, dy2/dt = y1  --> oscillating solution in y and w
        f4(t, y) = [-y[2], y[1]]
        sol4 = hcat(cos(t)-2*sin(t), 2*cos(t) + sin(t))

        context("The integrator initialises properly") do
            integrator = Ogre.GenericRK4
            s0 = integrator(f1, 0., t)
            s1 = integrator(f1, [0.], t)
            s2 = integrator(f4, [1., 2.], t)
            @fact start(s0) => (1, 0.0)
            @fact start(s1) => (1, [0.0])
            @fact start(s2) => (1, [1., 2.])
        end

        context("The dense solutions are as expected") do
            solver = Ogre.ode4_dense
            @fact solver(f1, 0., t) => roughly(sol1)
            @fact solver(f2, 0., t) => roughly(sol2)
            @fact solver(f3, 1., t) => roughly(sol3)
            @fact solver(f4, [1., 2.], t) => roughly(sol4, rtol = 1e-3)
        end
    end

    context("Planetary integrator test cases") do
        R_earth = Ogre.R_earth # m
        M_earth = Ogre.M_earth # kg
        surface_pressure = 1e5 # Pa

        NoTemp = Ogre.NoTemp
        WithTemp = Ogre.WithTemp

        dual_layer_eoses = [Ogre.mgsio3, Ogre.fe]
        tri_layer_eoses = [Ogre.h2o, Ogre.mgsio3, Ogre.fe]
        dual_layer_transitions = [0, 2/3, 1] .* M_earth
        tri_layer_transitions = [0, 1/6, 2/3, 1] * M_earth

        dual_layer_eos = Ogre.MassPiecewiseEOS(dual_layer_eoses,
                                              dual_layer_transitions)
        tri_layer_eos = Ogre.MassPiecewiseEOS(tri_layer_eoses,
                                              tri_layer_transitions)

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
                @fact (Ogre.R(M_earth, dual_layer_eos)
                       => roughly(R_earth, rtol=0.05))
                @fact (Ogre.R(M_earth, tri_layer_eos)
                       => roughly(R_earth, rtol=0.05))
            end
        end

        context("Solving for internal structure (tri-layer)") do
            actual_radius = Ogre.R(M_earth, tri_layer_eos)
            r_bracket = [actual_radius, actual_radius]
            bvs = Ogre.BoundaryValues(M_earth, actual_radius, surface_pressure)
            grid = linspace(M_earth, 0, Ogre.total_points)

            system = Ogre.PlanetSystem(M_earth, tri_layer_eos,
                                       bvs, grid, r_bracket)
            soln = Ogre.solve(system)

            context("Solution outputs are the correct type") do
                @fact isa(soln, Ogre.PlanetStructure{NoTemp}) => true
            end

            context("Solution boundaries match the boundary conditions") do
                @fact Ogre.mass(Ogre.centre(soln)) => 0
                @fact Ogre.mass(Ogre.surface(soln)) => M_earth
                @fact Ogre.pressure(Ogre.centre(soln)) =>
                      greater_than(Ogre.pressure(Ogre.surface(soln)))
                @fact Ogre.pressure(Ogre.surface(soln)) => surface_pressure
                @fact Ogre.radius(Ogre.centre(soln)) => less_than(100)
                @fact Ogre.radius(Ogre.centre(soln)) => greater_than(0)
                @fact Ogre.radius(Ogre.surface(soln)) => actual_radius
            end
        end

        context("Solving a single-layer solution two different ways") do
            f(P) = Ogre.mgsio3_func(P)
            f(P, T) = f(P)
            Cₚ = Ogre.HeatCapacity(T -> 1)
            eos_notemp = Ogre.SimpleEOS(NoTemp, f, "No temperature dependence")
            eos_withtemp = Ogre.SimpleEOS(WithTemp, f, "With temp dependence")

            actual_radius = Ogre.R(M_earth, eos_notemp)
            tempdep_radius = Ogre.R(M_earth, eos_withtemp, Cₚ)
            @fact tempdep_radius => actual_radius
        end

        context("Helper function tests") do
            system = Ogre.DefaultPlanetSystem(M_earth, dual_layer_eos)

            context("Determining when we're close to the centre of the planet") do
                @fact Ogre.hit_the_centre(-50) => true
                @fact Ogre.hit_the_centre(50) => false
                @fact Ogre.hit_the_centre(1e6) => false
                @fact Ogre.not_far_enough(1e6) => true
                @fact Ogre.not_far_enough(50) => false
            end
        end

        context("Temperature dependence") do
            R = R_earth
            R_bracket = [0, 2] * R_earth
            M = M_earth
            Psurf = 1e5 # Pa
            Tsurf = 300 # K
            h2o_heatcap_func(T::Real) = 4200 # J kg⁻¹ K⁻¹
            Cₚ = Ogre.HeatCapacity(h2o_heatcap_func)
            eos = Ogre.SimpleEOS(WithTemp, (P, T) -> 1000, "Test EOS")
            m_inner, m_outer = 0, M_earth
            solution_grid = linspace(m_outer, m_inner, 100)
            boundary = Ogre.BoundaryValues(M, R, Psurf, Tsurf)
            system = Ogre.PlanetSystem(M, eos, Cₚ, boundary, solution_grid,
                                       R_bracket)
            soln = Ogre.solve(system)

            #TODO: adjust R finding code to return the structure in parallel?

            mass = Ogre.mass
            radius = Ogre.radius
            pressure = Ogre.pressure
            temperature = Ogre.temperature
            centre = Ogre.centre
            surface = Ogre.surface

            @fact mass(centre(soln)) => 0
            @fact mass(surface(soln)) => M_earth
            @fact pressure(centre(soln)) => greater_than(pressure(surface(soln)))
            @fact pressure(surface(soln)) => surface_pressure
            @fact radius(surface(soln)) => R
            @fact temperature(centre(soln)) => greater_than(temperature(surface(soln)))
            @fact temperature(surface(soln)) => Tsurf
            @fact isa(soln, Ogre.PlanetStructure{WithTemp}) => true
        end
    end
end


