using ogre.common, ogre.eos
using FactCheck, PyCall

facts("Equation of state handling") do

    context("Piecewise equations of state") do
        transition_masses = [0, 1, 2, 3]
        eqn1(x) = x
        eqn2(x) = 2x
        eqn3(x) = 3x

        eos1 = SimpleEOS(eqn1, "eos1")
        eos2 = SimpleEOS(eqn2, "eos2")
        eos3 = SimpleEOS(eqn3, "eos3")

        piecewise = MassPiecewiseEOS([eos1, eos2, eos3], transition_masses)

        context("--get number of equations in each EOS") do
            @fact ogre.eos.n_eqs(eos1) => 1
            @fact ogre.eos.n_eqs(piecewise) => 3
        end

        context("--the PiecewiseEOS returns correct individual EOS") do
            @fact ogre.eos.get_layer_eos(piecewise, -1) => eos1
            @fact ogre.eos.get_layer_eos(piecewise, 0.5) => eos1
            @fact ogre.eos.get_layer_eos(piecewise, 1.5) => eos2
            @fact ogre.eos.get_layer_eos(piecewise, 2.5) => eos3
            @fact ogre.eos.get_layer_eos(piecewise, 5) => eos3
        end
    end

    context("TFD tests") do
        A = [55.845, 28.0855]
        Z = [26, 14]
        n = [1.6, 1]
        pressures = [0.1e6, 1e6, 10e6]                   # in Pa
        anticipated_densities = [5060, 7170, 13900]      # in kg/m3
        TFD_test(P) = ogre.eos.TFD(P, Z, A, n)
        @vectorize_1arg Real TFD_test

        @fact TFD_test(pressures) => roughly(anticipated_densities; rtol=0.05)
    end


end
