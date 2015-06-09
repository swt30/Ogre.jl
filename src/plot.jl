# PLOT.JL
# Plotting functions

export plot, phaseplot, plt, ticker

# Python plot setup
#------------------------------------------------------------------------------

using PyCall, PyPlot
@pyimport matplotlib.style as plotstyle
plotstyle.use("fivethirtyeight")
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.ticker as ticker

module conf
    "Common plotting configuration for structural plots"
    const kws = Dict(:linewidth=>2)
end

# Planet structure plotting
#------------------------------------------------------------------------------

import PyPlot.plot

"Label the mass coordinate x-axis"
function label_x_as_mass()
    xlabel(L"Mass enclosed / M$_\oplus$")
end

"Plot the temperature profile of a planet"
function plot_temperature_profile(soln::PlanetStructure{WithTemp}; kws...)
    x = vec(mass(soln)) / M_earth
    T = vec(temperature(soln))
    plot(x, T; conf.kws..., kws...)
    label_x_as_mass()
    ylabel("Temperature / K")
end

"Plot the radius profile of a planet"
function plot_radius_profile(soln::PlanetStructure; kws...)
    x = vec(mass(soln)) / M_earth
    r = vec(radius(soln)) / R_earth
    plot(x, r; conf.kws..., kws...)
    label_x_as_mass()
    ylabel(L"Radius / R$_\oplus$")
end

"Plot the pressure profile of a planet"
function plot_pressure_profile(soln::PlanetStructure; kws...)
    x = vec(mass(soln)) / M_earth
    P = vec(pressure(soln)) / 1e9
    plot(x, P; conf.kws..., kws...)
    label_x_as_mass()
    ylabel("Pressure / GPa")
end

"Plot the heat capacity profile of a planet"
function plot_heat_capacity_profile(soln::PlanetStructure{WithTemp}, 
    sys::PlanetSystem{WithTemp}; kws...)
    x = vec(mass(soln)) / M_earth
    P = vec(pressure(soln))
    T = vec(temperature(soln))
    thermal_gradient = sys.structure_equations[3]
    eos = thermal_gradient.eos
    heatcap = thermal_gradient.heatcap

    cₚ = map(heatcap, P, T)
    plot(x, cₚ; conf.kws..., kws...)
    label_x_as_mass()
    ylabel(L"Heat capacity, c$_p$ / J kg$^{-1}$ K$^{-1}$")
end

"Plot the thermal expansivity profile of a planet"
function plot_expansivity_profile(soln::PlanetStructure{WithTemp}, 
    sys::PlanetSystem{WithTemp}; kws...)
    x = vec(mass(soln)) / M_earth
    P = vec(pressure(soln))
    T = vec(temperature(soln))
    thermal_gradient = sys.structure_equations[3]
    eos = thermal_gradient.eos
    heatcap = thermal_gradient.heatcap

    αᵥ = map((P, T) -> thermal_expansivity(P, T, eos), P, T)
    plot(x, αᵥ; conf.kws..., kws...)
    label_x_as_mass()
    ylabel(L"Thermal expansivity, $\alpha$ / K$^{-1}$")
end

"Plot the density profile of a planet"
function plot_density_profile(soln::PlanetStructure{WithTemp}, 
    sys::PlanetSystem{WithTemp}; kws...)
    x = vec(mass(soln)) / M_earth
    P = vec(pressure(soln))
    T = vec(temperature(soln))
    thermal_gradient = sys.structure_equations[3]
    eos = thermal_gradient.eos
    heatcap = thermal_gradient.heatcap

    ρ = map(eos, P, T)
    plot(x, ρ; conf.kws..., kws...)
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

function plot(pb::PhaseBoundary; kws...)
    plot(pb.P, pb.T; kws...)
end

"Plot the P-T profile of a planet overlaid on the phase boundaries of water"
function phaseplot(soln::PlanetStructure{WithTemp}; kws...)
    P = vec(pressure(soln)) / 1e9
    T = vec(temperature(soln))

    figure()

    ax1 = subplot(111)
    xlabel("Pressure / GPa")
    ylabel("Temperature / K")
    xscale("log")
    yscale("log")

    plot_phases()
    plot(P, T, linewidth=2; kws...)

    xlim(xmin=1e-5, xmax=100)
    ylim(ymin=200, ymax=2000)

    tight_layout()
end

"Plot the phase boundaries of water in P-T space"
function plot_phases()
    for pb in Ogre.phase_boundaries
        new_pressure = pb.P/1e9
        adjusted_P = cpmod(pb, P=new_pressure)
        c = isa(pb, OtherPhaseBoundary) ? "Red" : "Black" 
        plot(adjusted_P, linewidth=1, color=c)
    end
end

