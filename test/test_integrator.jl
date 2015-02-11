import Ogre: Common, Integrator, Structure
using FactCheck

facts("Integrator tests") do

    context("Standard RK4 integrator test cases") do
        t = [0:0.1:1]

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
            integrator = Integrator.GenericRK4
            s0 = integrator(f1, 0., t)
            s1 = integrator(f1, [0.], t)
            s2 = integrator(f4, [1., 2.], t)
            @fact start(s0) => (1, 0.0)
            @fact start(s1) => (1, [0.0])
            @fact start(s2) => (1, [1., 2.])
        end

        context("The dense solutions are as expected") do
            solver = Integrator.ode4_dense
            @fact solver(f1, 0., t) => roughly(sol1)
            @fact solver(f2, 0., t) => roughly(sol2)
            @fact solver(f3, 1., t) => roughly(sol3)
            @fact solver(f4, [1., 2.], t) => roughly(sol4, rtol=1e-4)
        end
    end

    context("Planetary integrator test cases") do
        R_earth = 6.3781e6 # m
        M_earth = 5.972e24 # kg

        dual_layer_eoses = [Eos.mgsio3, Eos.fe]
        tri_layer_eoses = [Eos.h2o, Eos.mgsio3, Eos.fe]
        dual_layer_transitions = [0, 2/3, 1] .* M_earth
        tri_layer_transitions = [0, 1/6, 2/3, 1] * M_earth

        dual_layer_eos = Eos.MassPiecewiseEOS(dual_layer_eoses,
                                              dual_layer_transitions)
        tri_layer_eos = Eos.MassPiecewiseEOS(tri_layer_eoses,
                                              tri_layer_transitions)

        context("Solving for a radius") do
            context("We can match Madhu and Sara's results to within 1%") do
                test_masses = [0.5, 1.0, 5.0]
                fe = Eos.fe_seager
                mgsio3 = Eos.mgsio3_seager
                target_fe_radii = [0.628896, 0.76727, 1.17721]
                target_mgsio3_radii = [0.842621, 1.04239, 1.65182]

                fe_radii = Integrator.R(test_masses, fe, in_earth_units=true)
                mgsio3_radii = Integrator.R(test_masses, mgsio3, in_earth_units=true)

                @fact fe_radii => roughly(target_fe_radii, rtol=0.01)
                @fact mgsio3_radii => roughly(target_mgsio3_radii, rtol=0.01)
            end

            context("Dual- and tri-layer solutions work") do
                # accept a pretty large error since we're really just looking
                # to make sure that it produces the right type of solution
                @fact (Integrator.R(M_earth, dual_layer_eos)
                       => roughly(R_earth, rtol=0.1))
                @fact (Integrator.R(M_earth, tri_layer_eos)
                       => roughly(R_earth, rtol=0.1))
            end
        end

        context("Solving for internal structure (tri-layer)") do
            actual_radius = Integrator.R(M_earth, tri_layer_eos)
            R_bracket = [actual_radius, actual_radius]

            system = Integrator.setup_system(M_earth, tri_layer_eos, R_bracket)
            soln = Integrator.solve(system)

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
    end
end



