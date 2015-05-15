# PLOT.JL
# Plotting functions

export plot, plt, ticker

# Python plot setup
#------------------------------------------------------------------------------

using PyCall, PyPlot
@pyimport matplotlib.style as plotstyle
plotstyle.use("fivethirtyeight")
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.ticker as ticker

# Planet structure plotting
#------------------------------------------------------------------------------

import PyPlot.plot
function plot(soln::PlanetStructure{NoTemp})
    x = vec(mass(soln)) ./ M_earth
    r = vec(radius(soln)) ./ R_earth
    P = vec(pressure(soln)) ./ 1e9
    # the use of vec() is because these functions return array views
    # (SubArrays) which are apparently not popular with PyPlot

    ax1 = subplot(211)
    ylabel("Pressure / GPa")
    setp(ax1[:get_xticklabels](), visible=false)

    ax2 = subplot(212, sharex=ax1)
    xlabel(L"Mass enclosed / M$_\oplus$")
    ylabel(L"Radius / R$_\oplus$")

    for ax in gcf()[:get_axes]()
        ax[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 4))
    end

    ax1[:plot](x, P, linewidth=2)
    ax2[:plot](x, r, linewidth=2)
    tight_layout()
end

function plot(soln::PlanetStructure{WithTemp})
    x = vec(mass(soln)) ./ M_earth
    r = vec(radius(soln)) ./ R_earth
    P = vec(pressure(soln)) ./ 1e9
    T = vec(temperature(soln))

    ax1 = subplot(311)
    ylabel("Temperature / K")
    setp(ax1[:get_xticklabels](), visible=false)

    ax2 = subplot(312, sharex=ax1)
    ylabel("Pressure / GPa")
    setp(ax2[:get_xticklabels](), visible=false)

    ax3 = subplot(313, sharex=ax1)
    xlabel(L"Mass enclosed / M$_\oplus$")
    ylabel(L"Radius / R$_\oplus$")

    for ax in gcf()[:get_axes]()
        ax[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 4))
    end

    ax1[:plot](x, T, linewidth=2)
    ax2[:plot](x, P, linewidth=2)
    ax3[:plot](x, r, linewidth=2)
    tight_layout()
end


