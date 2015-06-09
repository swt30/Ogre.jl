module testresources

import Ogre

# ValueSets
value_set_no_temp = Ogre.ValueSet(1,2,3)
value_set_full = Ogre.ValueSet(1,2,3,4)

# Heat capacities
module heatcap
    import Ogre
    func = exp
    exponential = Ogre.HeatCapacity(Ogre.WithTemp, func)
    func2d(x, y) = exp(x) * exp(y)
    exponential2d = Ogre.HeatCapacity(Ogre.WithTempPressure, func2d)
    constant = Ogre.HeatCapacity(1)
end # module heatcap

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

module eos
    import Ogre

    # T inependent
    linear_f(x) = x
    squared_f(x) = x^2
    log_f(x) = log10(x)
    simple_fs = [linear_f, squared_f, log_f]
    simples = map(f -> Ogre.SimpleEOS(Ogre.NoTemp, f, ""), simple_fs)

    # piecewise
    piecewise = Ogre.MassPiecewiseEOS(simples, [0, 1, 2, 3])
    P_piecewise = Ogre.PressurePiecewiseEOS(simples, [0, 1, 2, 3])
    h2o_VII_seager_transition_pressures = [0, 44.3e9, 7686e9, 1e20]
    h2o_VII_seager_f(rho::Real) = Ogre.BME(rho, 1460., 23.7, 4.15) * 1e9
    h2o_seager_low = Ogre.InvPressureEOS(h2o_VII_seager_f,
                                         1e3, 1e8,
                                         "H2O (BME3) (Seager 2007)")
    # TODO: this line is fragile as it relies on the data directory - add data to test dir
    h2o_seager_dft = Ogre.load_interpolated_eos("$(Ogre.DATADIR)/tabulated/h2o-dft.eos")
    h2o_tfd_f(P::Real) = Ogre.TFD(P, [1, 8], [1.00794, 15.9994], [2., 1.])
    h2o_tfd = Ogre.SimpleEOS(Ogre.NoTemp, h2o_tfd_f, "H2O TFD")

    h2o_VII_seager_individual_eoses = [h2o_seager_low, h2o_seager_dft, h2o_tfd]
    h2o_VII_seager = Ogre.PressurePiecewiseEOS(h2o_VII_seager_individual_eoses,
                                               h2o_VII_seager_transition_pressures)

    # inverted
    make_inverted_eos(f::Function) = Ogre.InvPressureEOS(f, 1e5, 1e10,
                                                         "a test inverted EOS")
    simple_inverteds = map(make_inverted_eos, simple_fs)

    # T dependent
    PT_f(P, T) = P * T
    P2T2_f(P, T) = P^2 + T^2
    logPdivT_f(P, T) = log10(P / T)
    Tdep_fs= [PT_f, P2T2_f, logPdivT_f]
    Tdeps = map(f -> Ogre.SimpleEOS(Ogre.WithTemp, f, ""), Tdep_fs)

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

    interp_1D_lin = Ogre.lininterp(Plin, rholin_1D)
    interp_1D_log = Ogre.loginterp(Plog, rholog_1D)
    interp_2D_lin = Ogre.lininterp(Plin, Tlin, rholin_2D)
    interp_2D_log = Ogre.loginterp(Plog, Tlog, rholog_2D)
    interp_2D_lin_withnans = Ogre.lininterp(Plin, Tlin, rholin_2D_withnans; 
                                            suppress_warnings=true)
end # module eos

module planet
    module shared
        import Ogre
        const R_earth = Ogre.R_earth # m
        const M_earth = Ogre.M_earth # kg
        const atmospheric_pressure = 1e5 # Pa
        const surface_temperature = 300 # K
        const earth_mass_solution_grid = linspace(M_earth, 0, 100)
        const NoTemp = Ogre.NoTemp
        const WithTemp = Ogre.WithTemp

        export Ogre, R_earth, M_earth, atmospheric_pressure, surface_temperature,
               earth_mass_solution_grid, NoTemp, WithTemp
    end # module shared
    using .shared

    # temperature independent
    module tri_layer
        using ..shared
        eoses = [Ogre.h2o, Ogre.mgsio3, Ogre.fe]
        transitions = [0, 1/6, 2/3, 1] * M_earth
        eos = Ogre.MassPiecewiseEOS(eoses, transitions)
        actual_radius = Ogre.R(M_earth, eos)
        r_bracket = fill(actual_radius, 2)
        boundary_vals = Ogre.BoundaryValues(Ogre.M_earth, actual_radius, atmospheric_pressure)
        
        system = Ogre.PlanetSystem(M_earth, eos, boundary_vals, 
                               	   earth_mass_solution_grid, r_bracket)
    end # module tri_layer

    module dual_layer
        using ..shared
        eoses = [Ogre.mgsio3, Ogre.fe]
        transitions = [0, 2/3, 1] .* M_earth
        eos = Ogre.MassPiecewiseEOS(eoses, transitions)
    end # module dual_layer

    # temperature dependent
    module withtemp
        using ..shared
        R_bracket = [0, 2] * R_earth
        heatcap = Ogre.HeatCapacity(4200) # J K⁻¹ kg⁻¹
        eos = Ogre.SimpleEOS(WithTemp, (P, T) -> 1000 - 0.001T, "Test EOS")
        earth_surface_values = Ogre.BoundaryValues(M_earth, R_earth, 
                                                   atmospheric_pressure, surface_temperature)
        system = Ogre.PlanetSystem(M_earth, eos, heatcap, earth_surface_values, 
    						   		earth_mass_solution_grid, R_bracket)
    end # module withtemp
end # module planet

end # module testresources
