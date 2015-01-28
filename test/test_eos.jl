import Ogre: Common, Eos
using FactCheck

facts("Equation of state handling") do
    eqn1(x) = x
    eqn2(x) = x^2
    eqn3(x) = log10(x)
    simple_eos1 = Eos.SimpleEOS(eqn1, "a test EOS")
    simple_eos2 = Eos.SimpleEOS(eqn2, "a test EOS")
    simple_eos3 = Eos.SimpleEOS(eqn3, "a test EOS")

    context("Calling a simple EOS") do
        @fact Common.callfunc(simple_eos1, 42) => 42
    end
    @pending "Update this once call overrides are a thing" => ""

    context("Piecewise equations of state") do
        piecewise_eos = Eos.MassPiecewiseEOS([simple_eos1,
                                              simple_eos2,
                                              simple_eos3], [0, 1, 2, 3])

        context("get number of equations in each EOS") do
            @fact Eos.n_eqs(simple_eos1) => 1
            @fact Eos.n_eqs(piecewise_eos) => 3
        end

        context("the PiecewiseEOS returns correct individual EOS") do
            @fact Eos.get_layer_eos(piecewise_eos, -1)  => simple_eos1
            @fact Eos.get_layer_eos(piecewise_eos, 0.5) => simple_eos1
            @fact Eos.get_layer_eos(piecewise_eos, 1.5) => simple_eos2
            @fact Eos.get_layer_eos(piecewise_eos, 2.5) => simple_eos3
            @fact Eos.get_layer_eos(piecewise_eos, 5)   => simple_eos3
        end
    end

    context("TFD tests") do
        pressures = [0.1e6, 1e6, 10e6] # in bar
        pressures = pressures .* 1e5   # 1 bar = 1e5 Pa

        function test_TFD(A, Z, n, anticipated_densities)
            TFD(P) = Eos.TFD(P, Z, A, n)
            @vectorize_1arg Real TFD
            @fact TFD(pressures) => roughly(anticipated_densities,
                                                 rtol=0.01)
        end

        function test_TFD(A, Z, anticipated_densities)
            TFD(P) = Eos.TFD(P, Z, A)
            @vectorize_1arg Real TFD
            @fact TFD(pressures) => roughly(anticipated_densities,
                                                 rtol=0.01)
        end

        context("with a single element") do
            context("Fe") do
                test_TFD(55.845, 26, [5900, 8130, 15400])
            end
            context("Bi") do
                test_TFD(208.98, 83, [21800, 26200, 41000])
            end
        end

        context("with multiple elements") do
            context("TiO2") do
                test_TFD([47.867, 15.9994], [22, 8], [1., 2.],
                         [3190, 4920, 10400])
            end
            context("PbS") do
                test_TFD([207.2, 32.065], [82, 16], [12800, 17000, 29600])
            end
        end
    end

    context("EOS inversion") do
        make_inverted_eos(f::Function) = Eos.InvertedEOS(f, 1e5, 1e10,
                                                         "a test EOS")

        context("Known functions invert properly") do
            sample_density = 1e8
            testfuncs = [eqn1, eqn2, eqn3]
            eoses = map(make_inverted_eos, testfuncs)
            pressures = map(f -> f(sample_density), testfuncs)

            map(eoses, pressures) do eos, P
                @fact Common.callfunc(eos, P) => roughly(sample_density)
            end
            @pending "Update this once call overrides are a thing" => ""

        end

        context("Fail inverting outside the range specified") do
            f(rho) = rho
            eos = make_inverted_eos(f)
            sample_density = 1e12
            P = f(sample_density)

            @fact_throws Common.callfunc(f, P)
        end
    end

end
