module testresources

import Ogre
const DATADIR = Ogre.DATADIR

# Contains various resources to be used in testing

# ValueSets
value_set_no_temp = Ogre.ValueSet(1,2,3)
value_set_full = Ogre.ValueSet(1,2,3,4)

# EOS - T inependent
eos_linear_f(x) = x
eos_squared_f(x) = x^2
eos_log_f(x) = log10(x)
simple_eos_functions = [eos_linear_f, eos_squared_f, eos_log_f]
simple_eoses = map(f -> Ogre.SimpleEOS(Ogre.NoTemp, f, ""), simple_eos_functions)

# EOS - piecewise
simple_piecewise_EOS = Ogre.MassPiecewiseEOS(simple_eoses,
                                             [0, 1, 2, 3])
simple_P_piecewise_EOS = Ogre.PressurePiecewiseEOS(simple_eoses,
                                                   [0, 1, 2, 3])
h2o_VII_seager_transition_pressures = [0, 44.3e9, 7686e9, 1e20]
h2o_VII_seager_f(rho::Real) = Ogre.BME(rho, 1460., 23.7, 4.15) * 1e9
h2o_seager_low = Ogre.InvPressureEOS(h2o_VII_seager_f,
                                     1e3, 1e8,
                                     "H2O (BME3) (Seager 2007)")
# TODO: this line is fragile as it relies on the data directory - add data to test dir
h2o_seager_dft = Ogre.load_interpolated_eos("$DATADIR/tabulated/H2O (DFT).eos")
h2o_tfd_f(P::Real) = Ogre.TFD(P, [1, 8], [1.00794, 15.9994], [2., 1.])
h2o_tfd = Ogre.SimpleEOS(Ogre.NoTemp, h2o_tfd_f, "H2O TFD")

h2o_VII_seager_individual_eoses = [h2o_seager_low, h2o_seager_dft, h2o_tfd]
h2o_VII_seager = Ogre.PressurePiecewiseEOS(h2o_VII_seager_individual_eoses,
                                           h2o_VII_seager_transition_pressures)

# EOS - inverted
make_inverted_eos(f::Function) = Ogre.InvPressureEOS(f, 1e5, 1e10,
                                                     "a test inverted EOS")
simple_inv_eoses = map(make_inverted_eos, simple_eos_functions)

# EOS - temperature dependent
eos_PT_f(P, T) = P * T
eos_P2T2_f(P, T) = P^2 + T^2
eos_logPdivT_f(P, T) = log10(P / T)
simple_Tdep_eos_functions = [eos_PT_f, eos_P2T2_f, eos_logPdivT_f]
simple_Tdep_eoses = map(f -> Ogre.SimpleEOS(Ogre.WithTemp, f, ""), simple_Tdep_eos_functions)

end
