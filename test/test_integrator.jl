using Ogre: Common, Integrator
using FactCheck

facts("Integrator tests") do

    context("RK4 integrator produces correct results for basic cases") do
        t = [0:0.1:1]

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
        sol4 = [cos(t) - 2*sin(t), 2*cos(t) + sin(t)]

        solver = Integrator.ode4

        @pending solver(f1, 0., t) => roughly(f1)
        @pending solver(f2, 0., t) => roughly(f2)
        @pending solver(f3, 1., t) => roughly(f3)
        @pending solver(f4, [1., 2.], t) => roughly(f4)
    end

end

