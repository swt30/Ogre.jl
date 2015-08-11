include("header.jl")

facts("Equation of state (EOS) handling") do
    context("without temperature dependence") do
        eos_linear, eos_squared, eos_log = res.eos.simples

        context("Calling a simple EOS") do
            @fact eos_linear(42) --> 42
            @fact eos_squared(42) --> 42^2
        end

        context("Calling using ValueSets") do
            vs = res.value_set_no_temp
            @fact eos_linear(vs) --> 3
            @fact eos_squared(vs) --> 3^2
        end

        context("Piecewise EOS") do
            context("Get number of equations in each EOS") do
                @fact length(eos_linear) --> 1
                @fact length(res.eos.P_piecewise) --> 3
            end

            context("The PiecewiseEOS returns correct individual EOS") do
                peos = res.eos.piecewise
                @fact Ogre.get_layer_eos(peos, -1)  --> eos_linear
                @fact Ogre.get_layer_eos(peos, 0.5) --> eos_linear
                @fact Ogre.get_layer_eos(peos, 1.5) --> eos_squared
                @fact Ogre.get_layer_eos(peos, 2.5) --> eos_log
                @fact Ogre.get_layer_eos(peos, 5)   --> eos_log
            end

            context("Pressure-piecewise EOS calls match the individual EOS") do
                context("A really simple piecewise EOS") do
                    peos = res.eos.P_piecewise
                    @fact peos(0.5) --> eos_linear(0.5)
                    @fact peos(1.5) --> eos_squared(1.5)
                    @fact peos(2.5) --> eos_log(2.5)
                end

                context("A more complicated pressure-piecewise EOS") do
                    # We construct a pressure-piecewise EOS based on Sara
                    # Seager's water EOS in her 2007 paper, then attempt to draw
                    # from it to see if it matches our expectations

                    test_pressures = [1e8, 1e11, 1e14]
                    eoses = res.eos.h2o_VII_seager_individual_eoses
                    piecewise = res.eos.h2o_VII_seager

                    for (eos, P) in zip(eoses, test_pressures)
                        @fact piecewise(P) --> eos(P)
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
                @fact TFD(pressures) --> roughly(anticipated_densities, rtol=0.01)
            end

            function test_TFD(A, Z, anticipated_densities)
                TFD(P) = Ogre.TFD(P, Z, A)
                @vectorize_1arg Real TFD
                @fact TFD(pressures) --> roughly(anticipated_densities, rtol=0.01)
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
                    test_TFD([47.867, 15.9994], [22, 8], [1., 2.], [3190, 4920, 10400])
                end
                context("PbS") do
                    test_TFD([207.2, 32.065], [82, 16], [12800, 17000, 29600])
                end
            end
        end

        context("EOS inversion") do
            context("Known functions invert properly") do
                sample_density = 1e8
                testfuncs = res.eos.simple_fs
                eoses = res.eos.simple_inverteds
                pressures = map(f -> f(sample_density), testfuncs)

                map(eoses, pressures) do eos, P
                    @fact eos(P) --> roughly(sample_density)
                end
            end

            context("Fail inverting outside the range specified") do
                inv_linear = res.eos.simple_inverteds[1]
                sample_density = 1e12
                P = res.eos.linear_f(sample_density)

                @fact_throws inv_linear(P)
            end
        end

        context("Interpolated EOS") do
            context("match the actual function values") do
                context("for a simple EOS") do
                    lininterp = res.eos.interp_1D_lin
                    loginterp = res.eos.interp_1D_log
                    values_to_check = [1, 2, 4.7, 7.5, 8]
                    expected_results = [P^2 for P in values_to_check]

                    @fact loginterp(values_to_check) --> roughly(expected_results, rtol=0.01)
                    @fact lininterp(values_to_check) --> roughly(expected_results, rtol=0.01)
                end

                context("for a 2D EOS") do
                    lininterp2d = res.eos.interp_2D_lin
                    loginterp2d = res.eos.interp_2D_log

                    P_to_check = [1, 2, 3, 5.8, 9.1]
                    T_to_check = [10, 25, 40.3, 64.3, 94]
                    expected_results = [P^2 + T^2 for (P, T) in zip(P_to_check, T_to_check)]

                    @fact (lininterp2d(P_to_check, T_to_check) --> 
                           roughly(expected_results, rtol=0.01))
                    @fact (loginterp2d(P_to_check, T_to_check) --> 
                           roughly(expected_results, rtol=0.01))

                    context("NaN values error out") do
                        rnan = copy(res.eos.rholin_2D)
                        # set a value near the edge = NaN 
                        rnan[2] = NaN
                        @fact_throws Ogre.lininterp(res.eos.Plin, res.eos.Tlin, rnan)
                    end
                end
            end
        end
    end

    context("with temperature dependence") do
        teos = res.eos.Tdeps[1]
        context("Calling a simple EOS") do
            @fact teos(12, 3) --> 36
            @fact_throws teos[1](12)
        end

        context("Calling using ValueSets") do
            @fact teos(res.value_set_full) --> 12
            @fact_throws teos(res.value_set_no_temp)
        end
    end

    context("Phase boundaries") do
        context("have correct long-short code mappings") do
            iskeyof(dict) = k -> k in keys(dict)
            isvaluein(dict) = v -> v in values(dict)
            phasedict = Ogre.phase_mappings

            @fact "liquid" --> iskeyof(phasedict)
            @fact "L" --> isvaluein(phasedict)
            @fact "L" --> iskeyof(phasedict)
            @fact "ice VII" --> iskeyof(phasedict)
            @fact "VII" --> iskeyof(phasedict)
            @fact "VII" --> isvaluein(phasedict)
            @fact phasedict["liquid"] --> "L"
            @fact phasedict["ice VII"] --> "VII"
        end

        context("Have correct end values on the boundaries") do
            pb = Ogre.PhaseBoundary

            maxtemp(pb) = maximum(pb.T)
            mintemp(pb) = minimum(pb.T)
            maxpress(pb) = maximum(pb.P)
            minpress(pb) = minimum(pb.P)

            @fact_throws pb(:L, :II)
            @fact_throws pb("L", "II")

            @fact mintemp(pb(:L, :VII)) --> roughly(355.0, rtol=0.01)
            @fact maxpress(pb(:L, :VII)) --> roughly(400000e5, rtol=0.01)
            @fact mintemp(pb("L", :VII)) --> roughly(355.0, rtol=0.01)
            @fact minpress(pb("L", :VII)) --> roughly(22160e5, rtol=0.01)
            @fact mintemp(pb(:L, "VII")) --> roughly(355.0, rtol=0.01)
            @fact mintemp(pb("L", "VII")) --> roughly(355.0, rtol=0.01)
            @fact maxtemp(pb(:VII, :X)) --> roughly(1500, rtol=0.01)
        end
    end
end

