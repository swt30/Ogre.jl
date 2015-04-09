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
    x = vec(mass(soln))
    r = vec(radius(soln))
    P = vec(pressure(soln))
    # the use of vec() is because these functions return array views
    # (SubArrays) which are apparently not popular with PyPlot

    fig = figure()
    ax1 = subplot(211)
    ylabel("Pressure, P / Pa")
    setp(ax1[:get_xticklabels](), visible=false)

    ax2 = subplot(212, sharex=ax1)
    xlabel("Mass enclosed, m / kg")
    ylabel("Radius, r / m")

    for ax in fig[:get_axes]()
        ax[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 4))
    end

    ax1[:plot](x, P, linewidth=2)
    ax2[:plot](x, r, linewidth=2)
    tight_layout()
end

function plot(soln::PlanetStructure{WithTemp})
    x = vec(mass(soln))
    r = vec(radius(soln))
    P = vec(pressure(soln))
    T = vec(temperature(soln))

    fig = figure()
    ax1 = subplot(311)
    ylabel("Temperature, T / K")
    setp(ax1[:get_xticklabels](), visible=false)

    ax2 = subplot(312, sharex=ax1)
    ylabel("Pressure, P / Pa")
    setp(ax2[:get_xticklabels](), visible=false)

    ax3 = subplot(313, sharex=ax1)
    xlabel("Mass enclosed, m / kg")
    ylabel("Radius, r / m")

    for ax in fig[:get_axes]()
        ax[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 4))
    end

    ax1[:plot](x, T, linewidth=2)
    ax2[:plot](x, P, linewidth=2)
    ax3[:plot](x, r, linewidth=2)
    tight_layout()
end


