import Ogre: Common, Eos
using FactCheck

callfunc = Common.callfunc

facts("Equation of state (EOS) handling") do
    eqn1(x) = x
    eqn2(x) = x^2
    eqn3(x) = log10(x)
    simple_eos1 = Eos.SimpleEOS(eqn1, "a test EOS")
    simple_eos2 = Eos.SimpleEOS(eqn2, "a test EOS")
    simple_eos3 = Eos.SimpleEOS(eqn3, "a test EOS")

    context("Calling a simple EOS") do
        @fact callfunc(simple_eos1, 42) => 42
    end
    @pending "Update this once call overrides are a thing" => ""

    context("Piecewise EOS") do
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

        context("pressure-piecewise EOS calls match the individual EOS") do
            context("A really simple piecewise EOS") do
                P_piecewise_eos = Eos.PressurePiecewiseEOS([simple_eos1,
                                                            simple_eos2,
                                                            simple_eos3],
                                                            [0, 1, 2, 3])
                @fact callfunc(P_piecewise_eos, 0.5) => callfunc(simple_eos1, 0.5)
                @fact callfunc(P_piecewise_eos, 1.5) => callfunc(simple_eos2, 1.5)
                @fact callfunc(P_piecewise_eos, 2.5) => callfunc(simple_eos3, 2.5)
            end

            context("A more complicated pressure-piecewise EOS") do
                # We construct a pressure-piecewise EOS based on Sara Seager's
                # water EOS in her 2007 paper, then attempt to draw from it to see
                # if it matches our expectations
                transition_pressures = [0, 44.3e9, 7686e9, 1e20]

                h2o_VII_seager_func(rho::Real) = Eos.BME(rho, 1460., 23.7, 4.15) * 1e9
                h2o_seager_low = Eos.InvertedEOS(h2o_VII_seager_func, 1e3, 1e8,
                                                 "H2O (BME3) (Seager 2007)")
                h2o_seager_dft = Eos.load_interpolated_eos("data/tabulated/H2O (DFT).eos")
                h2o_tfd_func(P::Real) = Eos.TFD(P, [1, 8], [1.00794, 15.9994], [2., 1.])
                h2o_tfd = Eos.SimpleEOS(h2o_tfd_func, "H2O TFD")

                eoses = [h2o_seager_low, h2o_seager_dft, h2o_tfd]
                P_piecewise_eos = Eos.PressurePiecewiseEOS(eoses,
                                                           transition_pressures)

                test_pressures = [1e8, 1e11, 1e14]

                for (eos, P) in zip(eoses, test_pressures)
                    @fact callfunc(P_piecewise_eos, P) => callfunc(eos, P)
                end
            end
            @pending "Update this once call overrides are a thing" => ""
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
                @fact callfunc(eos, P) => roughly(sample_density)
            end
            @pending "Update this once call overrides are a thing" => ""

        end

        context("Fail inverting outside the range specified") do
            f(rho) = rho
            eos = make_inverted_eos(f)
            sample_density = 1e12
            P = f(sample_density)

            @fact_throws callfunc(f, P)
            @pending "Update this once call overrides are a thing" => ""
        end
    end

    context("Interpolated EOS") do
        context("match the actual function values") do
            context("for a simple EOS") do
                Plin = linspace(1, 10)
                Plog = logspace(0, 1)
                rholin = Plin.^2
                rholog = Plog.^2
                log_interpolation = Eos.loginterp(Plog, rholog)
                lin_interpolation = Eos.lininterp(Plin, rholin)

                @fact log_interpolation(1) => roughly(1^2)
                @fact log_interpolation(2) => roughly(2^2)
                @fact log_interpolation(4.7) => roughly(4.7^2)
                @fact log_interpolation(7.5) => roughly(7.5^2)
                @fact log_interpolation(8) => roughly(8^2)

                @fact lin_interpolation(1) => roughly(1^2)
                @fact lin_interpolation(2) => roughly(2^2)
                @fact lin_interpolation(4.7) => roughly(4.7^2)
                @fact lin_interpolation(7.5) => roughly(7.5^2)
                @fact lin_interpolation(10) => roughly(10^2)
            end
        end
    end

end
