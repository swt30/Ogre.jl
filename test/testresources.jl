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
h2o_seager_dft = Ogre.load_interpolated_eos("$DATADIR/tabulated/h2o-dft.eos")
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

# EOS - interpolated
Plin = linspace(1,10)
Plog = logspace(0,1)
Tlin = linspace(10,100)
Tlog = logspace(1,2)
rholin_1D = map(P->P^2, Plin)
rholog_1D = map(P->P^2, Plog)
rholin_2D = Float64[P^2 + T^2 for P in Plin, T in Tlin]
rholog_2D = Float64[P^2 + T^2 for P in Plog, T in Tlog]
rholin_2D_withnans = copy(rholin_2D)
rholin_2D_withnans[2] = NaN

eos_interpolated_1D_linear = Ogre.lininterp(Plin, rholin_1D)
eos_interpolated_1D_log = Ogre.loginterp(Plog, rholog_1D)
eos_interpolated_2D_linear = Ogre.lininterp(Plin, Tlin, rholin_2D)
eos_interpolated_2D_log = Ogre.loginterp(Plog, Tlog, rholog_2D)
eos_interpolated_2D_linear_withnans = Ogre.lininterp(Plin, Tlin, rholin_2D_withnans; 
                                                     suppress_warnings=true)

# Heat capacities
heat_capacity_function = exp
exponential_heat_capacity = Ogre.HeatCapacity(Ogre.WithTemp, heat_capacity_function)

# Simple ODEs
# dy/dt = 6  -->  y = 6t (y0=0)
ode_1(t, y) = 6
ode_1_sol(t) = 6.*t
# dy/dt = 2t -->  t = t^2 + y0 (y0=0)
ode_2(t, y) = 2t
ode_2_sol(t) = t.^2
# dy/dt = y  -->  y = y0 e^t (y0=1)
ode_3(t, y) = y
ode_3_sol(t) = exp(t)
# dy1/dt = -y2, dy2/dt = y1  --> oscillating solution in y and w
ode_4(t, y) = [-y[2], y[1]]
ode_4_sol(t) = hcat(cos(t)-2*sin(t), 2*cos(t) + sin(t))

simple_odes = [ode_1, ode_2, ode_3, ode_4]
simple_ode_solutions = [ode_1_sol, ode_2_sol, ode_3_sol, ode_4_sol]

# Planetary structures
const R_earth = Ogre.R_earth # m
const M_earth = Ogre.M_earth # kg
const atmospheric_pressure = 1e5 # Pa
const surface_temperature = 300 # K

NoTemp = Ogre.NoTemp
WithTemp = Ogre.WithTemp

# temperature independent
dual_layer_eoses = [Ogre.mgsio3, Ogre.fe]
tri_layer_eoses = [Ogre.h2o, Ogre.mgsio3, Ogre.fe]
dual_layer_transitions = [0, 2/3, 1] .* M_earth
tri_layer_transitions = [0, 1/6, 2/3, 1] * M_earth
dual_layer_eos = Ogre.MassPiecewiseEOS(dual_layer_eoses,
                                      dual_layer_transitions)
tri_layer_eos = Ogre.MassPiecewiseEOS(tri_layer_eoses,
                                      tri_layer_transitions)

tri_layer_planet_actual_radius = Ogre.R(M_earth, tri_layer_eos)
tri_layer_planet_r_bracket = fill(tri_layer_planet_actual_radius, 2)
tri_layer_planet_boundary_vals = Ogre.BoundaryValues(Ogre.M_earth, 
								                     tri_layer_planet_actual_radius, 
								                     atmospheric_pressure)
earth_mass_planet_solution_grid = linspace(Ogre.M_earth, 0)
tri_layer_planet = Ogre.PlanetSystem(Ogre.M_earth, tri_layer_eos,
                           			 tri_layer_planet_boundary_vals, 
                           			 earth_mass_planet_solution_grid, 
                           			 tri_layer_planet_r_bracket)

# temperature dependent
R_bracket = [0, 2] * R_earth
simple_h2o_heatcap_func(T::Real) = 4200 # J kg⁻¹ K⁻¹
simple_h2o_heatcap = Ogre.HeatCapacity(simple_h2o_heatcap_func)
simple_constant_eos = Ogre.SimpleEOS(Ogre.WithTemp, (P, T) -> 1000, "Test EOS")
earth_surface_values = Ogre.BoundaryValues(M_earth, R_earth, atmospheric_pressure, 300)
simple_temperature_dependent_system = Ogre.PlanetSystem(M_earth, simple_constant_eos, 
														simple_h2o_heatcap, 
						   								earth_surface_values, 
						   								earth_mass_planet_solution_grid,
						   								R_bracket)

end
