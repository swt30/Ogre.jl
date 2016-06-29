# Plotting functions

using PyCall, PyPlot
import PyPlot: plot


# Plotting planetary structures

"Label the mass coordinate x-axis"
function label_x_as_mass()
    xlabel(L"Mass enclosed / M$_\oplus$")
end

"Plot the temperature profile of a planet"
function plot_temperature_profile(soln::PlanetStructure{WithTemp}; kws...)
    M = vec(mass(soln))
    T = vec(temperature(soln))
    plot(M/M_earth, T/K; kws...)
    label_x_as_mass()
    ylabel("Temperature / K")
end

"Plot the radius profile of a planet"
function plot_radius_profile(soln::PlanetStructure; kws...)
    M = vec(mass(soln))
    r = vec(radius(soln))
    plot(M/M_earth, r/R_earth; kws...)
    label_x_as_mass()
    ylabel(L"Radius / R$_\oplus$")
end

"Plot the pressure profile of a planet"
function plot_pressure_profile(soln::PlanetStructure; kws...)
    M = vec(mass(soln))
    P = vec(pressure(soln))
    plot(M/M_earth, P/GPa; kws...)
    label_x_as_mass()
    ylabel("Pressure / GPa")
end

"Plot the heat capacity profile of a planet"
function plot_heat_capacity_profile(soln::PlanetStructure{WithTemp},
                                    sys::PlanetSystem{WithTemp}; kws...)
    M = vec(mass(soln))
    P = vec(pressure(soln))
    T = vec(temperature(soln))
    thermal_gradient = sys.structure_equations[3]
    eos = thermal_gradient.eos
    heatcap = thermal_gradient.heatcap

    cₚ = map(heatcap, P, T)
    plot(M/M_earth, cₚ/(J/kg/K); kws...)
    label_x_as_mass()
    ylabel(L"Heat capacity, c$_p$ / J kg$^{-1}$ K$^{-1}$")
end

"Plot the thermal expansivity profile of a planet"
function plot_expansivity_profile(soln::PlanetStructure{WithTemp},
                                  sys::PlanetSystem{WithTemp}; kws...)
    M = vec(mass(soln))
    P = vec(pressure(soln))
    T = vec(temperature(soln))
    thermal_gradient = sys.structure_equations[3]
    eos = thermal_gradient.eos
    heatcap = thermal_gradient.heatcap

    αᵥ = map(M, P, T) do M, P, T
        vs = ValueSet(M, NaN, P, T)
        Ogre.thermexp(eos, vs)
    end

    plot(M/M_earth, αᵥ/(1/K); kws...)
    label_x_as_mass()
    ylabel(L"Thermal expansivity, $\alpha$ / K$^{-1}$")
end

"Plot the density profile of a planet"
function plot_density_profile(soln::PlanetStructure{WithTemp},
                              sys::PlanetSystem{WithTemp}; kws...)
    M = vec(mass(soln))
    P = vec(pressure(soln))
    T = vec(temperature(soln))
    thermal_gradient = sys.structure_equations[3]
    eos = thermal_gradient.eos
    heatcap = thermal_gradient.heatcap

    ρ = map((M, P, T) -> eos(Ogre.ValueSet(M, NaN, P, T)), M, P, T)
    plot(M/M_earth, ρ/(kg/m^3); kws...)
    label_x_as_mass()
    ylabel(L"Density, $\rho$ / kg m$^{-3}$")
end

function plot(soln::PlanetStructure{NoTemp}; kws...)
    fig, axes = subplots(2, 1, sharex=true)

    sca(axes[1])
    plot_pressure_profile(soln)
    setp(axes[1][:get_xticklabels](), visible=false)
    xlabel("")

    sca(axes[2])
    plot_radius_profile(soln)

    for ax in axes
        ax[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 4))
    end

    tight_layout()
end

function plot(soln::PlanetStructure{WithTemp}; kws...)
    fig, axes = subplots(3, 1, sharex=true)

    sca(axes[1])
    plot_temperature_profile(soln)

    sca(axes[2])
    plot_pressure_profile(soln)

    sca(axes[3])
    plot_radius_profile(soln)

    for ax in axes[1:2]
        setp(ax[:get_xticklabels](), visible=false)
        ax[:set_xlabel]("")
    end

    for ax in axes
        ax[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 4))
    end

    tight_layout()
end

function plot(soln::PlanetStructure{WithTemp}, sys::PlanetSystem{WithTemp}; kws...)
    fig, axes = subplots(3, 2, sharex=true, figsize=(14,9))

    sca(axes[1])
    plot_temperature_profile(soln)

    sca(axes[2])
    plot_pressure_profile(soln)

    sca(axes[3])
    plot_radius_profile(soln)

    sca(axes[4])
    plot_density_profile(soln, sys)

    sca(axes[5])
    plot_expansivity_profile(soln, sys)

    sca(axes[6])
    plot_heat_capacity_profile(soln, sys)

    for ax in axes[[1,2,4,5]]
        setp(ax[:get_xticklabels](), visible=false)
        ax[:set_xlabel]("")
    end

    for ax in axes
        ax[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 4))
    end
end

function plot(pb::PhaseBoundary; kwargs...)
    plot(pb.P/Pa, pb.T/K; kwargs...)
end

"Plot the P-T profile of a planet overlaid on the phase boundaries of water"
function phaseplot(soln::PlanetStructure{WithTemp}; kwargs...)
    P = vec(pressure(soln))
    T = vec(temperature(soln))

    figure()

    ax1 = subplot(111)
    xlabel("Pressure / Pa")
    ylabel("Temperature / K")
    xscale("log")
    yscale("log")

    plot_phases()
    plot(P/Pa, T/K, kwargs...)

    xlim(xmin=1e4, xmax=1e12)
    ylim(ymin=200, ymax=10000)

    tight_layout()
end

"Plot the phase boundaries of water in P-T space"
function plot_phases()
    boundaries = WaterData.load_phase_boundaries()["boundaries"]
    iapws = boundaries["iapws"]
    plot(iapws.P/Pa, iapws.T/K, linewidth=1, color="Red")
    for pb in boundaries["dunaeva"]
        plot(pb.P/Pa, pb.T/K, linewidth=1, color="Black")
    end
end
