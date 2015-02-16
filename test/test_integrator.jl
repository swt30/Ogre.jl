import Ogre
using FactCheck

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
        R_earth = 6.3781e6 # m
        M_earth = 5.972e24 # kg

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
                       => roughly(R_earth, rtol=0.1))
                @fact (Ogre.R(M_earth, tri_layer_eos)
                       => roughly(R_earth, rtol=0.1))
            end
        end

        context("Solving for internal structure (tri-layer)") do
            actual_radius = Ogre.R(M_earth, tri_layer_eos)
            R_bracket = [actual_radius, actual_radius]

            system = Ogre.setup_planet(M_earth, tri_layer_eos, R_bracket)
            soln = Ogre.solve(system)

            context("Solution outputs are the correct shape") do
                @fact ndims(soln.m) => 1
                @fact ndims(soln.y) => 2
            end

            context("Solution boundaries match the boundary conditions") do
                surface_pressure = 1e5
                r, P = soln.y[:, 1] , soln.y[:, 2]
                @fact soln.m[1] => M_earth
                @fact soln.m[end] => 0
                @fact P[1] => surface_pressure
                @fact r[1] => actual_radius
                @fact r[end] => less_than(100)
                @fact r[end] => greater_than(0)
            end
        end

        context("Helper function tests") do
            R_bracket = [0, 15] * R_earth
            system = Ogre.setup_planet(M_earth, dual_layer_eos, R_bracket)

            context("Determining when we're close to the centre of the planet") do
                @fact Ogre.hit_the_centre(-50) => true
                @fact Ogre.hit_the_centre(50) => false
                @fact Ogre.hit_the_centre(1e6) => false
                @fact Ogre.not_far_enough(1e6) => true
                @fact Ogre.not_far_enough(50) => false
            end
        end
    end
end



