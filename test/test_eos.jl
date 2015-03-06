include("header.jl")

# simplify a couple of common variables
const DATADIR = Ogre.DATADIR

facts("Equation of state (EOS) handling") do
    context("without temperature dependence") do
        # general setup of some simple EOSes
        eqn1(x) = x
        eqn2(x) = x^2
        eqn3(x) = log10(x)
        simple_eos1 = Ogre.SimpleEOS(Ogre.notemp, eqn1, "a test EOS")
        simple_eos2 = Ogre.SimpleEOS(Ogre.notemp, eqn2, "a test EOS")
        simple_eos3 = Ogre.SimpleEOS(Ogre.notemp, eqn3, "a test EOS")

        context("Calling a simple EOS") do
            @fact simple_eos1(42) => 42
        end

        context("Calling using ValueSets") do
            vs = Ogre.ValueSet(1,2,3) # mass, radius, pressure
            @fact simple_eos1(vs) => eqn1(3)
        end

        context("Piecewise EOS") do
            # stitch together those simple EOSes
            piecewise_eos = Ogre.MassPiecewiseEOS([simple_eos1,
                                                  simple_eos2,
                                                  simple_eos3], [0, 1, 2, 3])

            context("Get number of equations in each EOS") do
                @fact length(simple_eos1) => 1
                @fact length(piecewise_eos) => 3
            end

            context("The PiecewiseEOS returns correct individual EOS") do
                @fact Ogre.get_layer_eos(piecewise_eos, -1)  => simple_eos1
                @fact Ogre.get_layer_eos(piecewise_eos, 0.5) => simple_eos1
                @fact Ogre.get_layer_eos(piecewise_eos, 1.5) => simple_eos2
                @fact Ogre.get_layer_eos(piecewise_eos, 2.5) => simple_eos3
                @fact Ogre.get_layer_eos(piecewise_eos, 5)   => simple_eos3
            end

            context("Pressure-piecewise EOS calls match the individual EOS") do
                context("A really simple piecewise EOS") do
                    P_piecewise_eos = Ogre.PressurePiecewiseEOS([simple_eos1,
                                                                simple_eos2,
                                                                simple_eos3],
                                                                [0, 1, 2, 3])
                    @fact P_piecewise_eos(0.5) => simple_eos1(0.5)
                    @fact P_piecewise_eos(1.5) => simple_eos2(1.5)
                    @fact P_piecewise_eos(2.5) => simple_eos3(2.5)
                end

                context("A more complicated pressure-piecewise EOS") do
                    # We construct a pressure-piecewise EOS based on Sara Seager's
                    # water EOS in her 2007 paper, then attempt to draw from it to see
                    # if it matches our expectations

                    transition_pressures = [0, 44.3e9, 7686e9, 1e20]

                    h2o_VII_seager_func(rho::Real) = Ogre.BME(rho, 1460., 23.7, 4.15) * 1e9
                    h2o_seager_low = Ogre.InvPressureEOS(h2o_VII_seager_func,
                                                         1e3, 1e8,
                                                         "H2O (BME3) (Seager 2007)")

                    # TODO: this line is fragile as it relies on the data directory - add data to test dir
                    h2o_seager_dft = Ogre.load_interpolated_eos("$DATADIR/tabulated/H2O (DFT).eos")
                    h2o_tfd_func(P::Real) = Ogre.TFD(P, [1, 8], [1.00794, 15.9994], [2., 1.])
                    h2o_tfd = Ogre.SimpleEOS(Ogre.notemp, h2o_tfd_func, "H2O TFD")

                    eoses = [h2o_seager_low, h2o_seager_dft, h2o_tfd]
                    P_piecewise_eos = Ogre.PressurePiecewiseEOS(eoses,
                                                               transition_pressures)

                    test_pressures = [1e8, 1e11, 1e14]

                    for (eos, P) in zip(eoses, test_pressures)
                        @fact P_piecewise_eos(P) => eos(P)
                    end
                end
            end
        end

        context("TFD tests") do
            # our sample pressures
            pressures = [0.1e6, 1e6, 10e6] # in bar
            pressures .*= 1e5              # 1 bar = 1e5 Pa

            function test_TFD(A, Z, n, anticipated_densities)
                TFD(P) = Ogre.TFD(P, Z, A, n)
                @vectorize_1arg Real TFD
                @fact TFD(pressures) => roughly(anticipated_densities,
                                                     rtol=0.01)
            end

            function test_TFD(A, Z, anticipated_densities)
                TFD(P) = Ogre.TFD(P, Z, A)
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
            make_inverted_eos(f::Function) = Ogre.InvPressureEOS(f, 1e5, 1e10,
                                                                 "a test EOS")

            context("Known functions invert properly") do
                sample_density = 1e8
                testfuncs = [eqn1, eqn2, eqn3]
                eoses = map(make_inverted_eos, testfuncs)
                pressures = map(f -> f(sample_density), testfuncs)

                map(eoses, pressures) do eos, P
                    @fact eos(P) => roughly(sample_density)
                end
            end

            context("Fail inverting outside the range specified") do
                f(rho) = rho
                eos = make_inverted_eos(f)
                sample_density = 1e12
                P = f(sample_density)

                @fact_throws eos(P)
            end
        end

        context("Interpolated EOS") do
            context("match the actual function values") do
                context("for a simple EOS") do
                    Plin = linspace(1, 10)
                    Plog = logspace(0, 1)
                    rholin = Plin.^2
                    rholog = Plog.^2
                    log_interpolation = Ogre.loginterp(Plog, rholog)
                    lin_interpolation = Ogre.lininterp(Plin, rholin)

                    values_to_check = [1, 2, 4.7, 7.5, 8]
                    expected_results = values_to_check.^2

                    function check_eos{T<:Real}(eos::Function,
                                                inputs::Vector{T},
                                                outputs::Vector{T})
                        @assert length(inputs) == length(outputs)
                        map(inputs, outputs) do input, output
                            @fact eos(input) => roughly(output)
                        end
                    end

                    check_eos(log_interpolation, values_to_check, expected_results)
                    check_eos(lin_interpolation, values_to_check, expected_results)
                end
            end
        end
    end

    context("with temperature dependence") do
        # general setup of some simple EOSes
        eqn1(P, T) = P * T
        eqn2(P, T) = P^2 + T^2
        eqn3(P, T) = log10(P / T)
        simple_eos1 = Ogre.PressureTempEOS(eqn1, "a test EOS")
        simple_eos2 = Ogre.PressureTempEOS(eqn2, "a test EOS")
        simple_eos3 = Ogre.PressureTempEOS(eqn3, "a test EOS")

        context("Calling a simple EOS") do
            @fact simple_eos1(12, 3) => 36
            @fact_throws simple_eos1(12)
        end

        context("Calling using ValueSets") do
            vs = Ogre.ValueSet(1,2,3,4) # mass, radius, pressure, temperature
            vs2 = Ogre.ValueSet(1,2,3)
            @fact simple_eos1(vs) => eqn1(3,4)
            @fact_throws simple_eos1(vs2)
        end
    end
end

