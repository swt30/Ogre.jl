using ogre.common
using ogre.eos
using FactCheck

context("Piecewise equation of state handling") do
    transition_masses = [0, 1, 2, 3]
    eqn1(x) = x
    eqn2(x) = 2x
    eqn3(x) = 3x

    eos1 = SimpleEOS(eqn1)
    eos2 = SimpleEOS(eqn2)
    eos3 = SimpleEOS(eqn3)

    piecewise = PiecewiseEOS([eos1, eos2, eos3], transition_masses)

    # number of equations in the EOS
    @fact ogre.eos.n_eqs(eos1) => 1
    @fact ogre.eos.n_eqs(piecewise) => 3

    # the PiecewiseEOS returns correct individual EOS
    @fact ogre.eos.get_layer_eos(piecewise, -1) => eos1
    @fact ogre.eos.get_layer_eos(piecewise, 0.5) => eos1
    @fact ogre.eos.get_layer_eos(piecewise, 1.5) => eos2
    @fact ogre.eos.get_layer_eos(piecewise, 2.5) => eos3
    @fact ogre.eos.get_layer_eos(piecewise, 5) => eos3
end
