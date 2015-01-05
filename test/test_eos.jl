using Ogre: Common, Eos
using FactCheck

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
            @fact Ogre.Eos.n_eqs(eos1) => 1
            @fact Ogre.Eos.n_eqs(piecewise) => 3
        end

        context("--the PiecewiseEOS returns correct individual EOS") do
            @fact Ogre.Eos.get_layer_eos(piecewise, -1) => eos1
            @fact Ogre.Eos.get_layer_eos(piecewise, 0.5) => eos1
            @fact Ogre.Eos.get_layer_eos(piecewise, 1.5) => eos2
            @fact Ogre.Eos.get_layer_eos(piecewise, 2.5) => eos3
            @fact Ogre.Eos.get_layer_eos(piecewise, 5) => eos3
        end
    end

    context("TFD tests") do

        pressures = [0.1e6, 1e6, 10e6] # in bar
        pressures = pressures .* 1e5   # 1 bar = 1e5 Pa

        context("--with a single element") do
            A = 56.
            Z = 26

            anticipated_densities = [5900, 8130, 15400]
            TFD_test(P) = Ogre.Eos.TFD(P, Z, A)
            @vectorize_1arg Real TFD_test

            @fact TFD_test(pressures) => roughly(anticipated_densities,
                                                 rtol=0.02)
        end

        context("--with multiple elements") do
            A = [56., 28.]
            Z = [26, 14]
            n = [2.1, 1]

            anticipated_densities = [5060, 7170, 13900]      # in kg/m3
            TFD_test(P) = Ogre.Eos.TFD(P, Z, A, n)
            @vectorize_1arg Real TFD_test

            @fact TFD_test(pressures) => roughly(anticipated_densities,
                                                 rtol=0.02)
        end
    end


end
