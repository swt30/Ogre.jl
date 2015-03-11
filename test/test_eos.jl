include("header.jl")

facts("Equation of state (EOS) handling") do
    context("without temperature dependence") do
        eos_linear, eos_squared, eos_log = res.simple_eoses

        context("Calling a simple EOS") do
            @fact eos_linear(42) => 42
            @fact eos_squared(42) => 42^2
        end

        context("Calling using ValueSets") do
            vs = res.value_set_no_temp
            @fact eos_linear(vs) => 3
            @fact eos_squared(vs) => 9
        end

        context("Piecewise EOS") do
            context("Get number of equations in each EOS") do
                @fact length(eos_linear) => 1
                @fact length(res.simple_piecewise_EOS) => 3
            end

            context("The PiecewiseEOS returns correct individual EOS") do
                peos = res.simple_piecewise_EOS
                @fact Ogre.get_layer_eos(peos, -1)  => eos_linear
                @fact Ogre.get_layer_eos(peos, 0.5) => eos_linear
                @fact Ogre.get_layer_eos(peos, 1.5) => eos_squared
                @fact Ogre.get_layer_eos(peos, 2.5) => eos_log
                @fact Ogre.get_layer_eos(peos, 5)   => eos_log
            end

            context("Pressure-piecewise EOS calls match the individual EOS") do
                context("A really simple piecewise EOS") do
                    peos = res.simple_P_piecewise_EOS
                    @fact peos(0.5) => eos_linear(0.5)
                    @fact peos(1.5) => eos_squared(1.5)
                    @fact peos(2.5) => eos_log(2.5)
                end

                context("A more complicated pressure-piecewise EOS") do
                    # We construct a pressure-piecewise EOS based on Sara Seager's
                    # water EOS in her 2007 paper, then attempt to draw from it to see
                    # if it matches our expectations

                    test_pressures = [1e8, 1e11, 1e14]
                    eoses = res.h2o_VII_seager_individual_eoses
                    piecewise = res.h2o_VII_seager

                    for (eos, P) in zip(eoses, test_pressures)
                        @fact piecewise(P) => eos(P)
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
            context("Known functions invert properly") do
                sample_density = 1e8
                testfuncs = res.simple_eos_functions
                eoses = res.simple_inv_eoses
                pressures = map(f -> f(sample_density), testfuncs)

                map(eoses, pressures) do eos, P
                    @fact eos(P) => roughly(sample_density)
                end
            end

            context("Fail inverting outside the range specified") do
                inv_linear = res.simple_inv_eoses[1]
                sample_density = 1e12
                P = res.eos_linear_f(sample_density)

                @fact_throws inv_linear(P)
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
        teos = res.simple_Tdep_eoses[1]
        context("Calling a simple EOS") do
            @fact teos(12, 3) => 36
            @fact_throws teos[1](12)
        end

        context("Calling using ValueSets") do
            @fact teos(res.value_set_full) => 12
            @fact_throws teos(res.value_set_no_temp)
        end
    end
end

