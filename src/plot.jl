# PLOT.JL
# Plotting functions

export plot, plt, ticker

# Python plot setup
#------------------------------------------------------------------------------

using PyCall
style = pyimport("matplotlib.style")
plt = pyimport("matplotlib.pyplot")
ticker = pyimport("matplotlib.ticker")
style[:use]("fivethirtyeight") # clean up this once dot-overloading is allowed

# export the generic plot function for use in interactive work
function plot(args...; kwargs...)
    plt[:plot](args...; kwargs...)
end

# Planet structure plotting
#------------------------------------------------------------------------------

function plot(soln::PlanetStructure{NoTemp})
    x = vec(mass(soln))
    r = vec(radius(soln))
    P = vec(pressure(soln))
    # the use of vec() is because these functions return array slices
    # (SubArrays) which are apparently not popular with PyPlot

    fig = plt[:figure]()
    ax = fig[:add_subplot](211)
    ax[:set_ylabel]("Pressure, P / Pa")
    plt[:setp](ax[:get_xticklabels](), visible=false)

    ax2 = fig[:add_subplot](212, sharex=ax)
    ax2[:set_xlabel]("Mass enclosed, m / kg")
    ax2[:set_ylabel]("Radius, r / m")

    for axes in fig[:get_axes]()
        axes[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 1))
    end

    fig[:tight_layout]()

    ax[:plot](x, P, linewidth=2)
    ax2[:plot](x, r, linewidth=2)

    plt[:show]()
end

function plot(soln::PlanetStructure{WithTemp})
    x = vec(mass(soln))
    r = vec(radius(soln))
    P = vec(pressure(soln))
    T = vec(temperature(soln))

    fig = plt[:figure]()
    ax = fig[:add_subplot](311)
    ax[:set_ylabel]("Temperature, T / K")
    plt[:setp](ax[:get_xticklabels](), visible=false)

    ax2 = fig[:add_subplot](312, sharex=ax)
    ax2[:set_ylabel]("Pressure, P / Pa")
    plt[:setp](ax2[:get_xticklabels](), visible=false)

    ax3 = fig[:add_subplot](313, sharex=ax)
    ax3[:set_xlabel]("Mass enclosed, m / kg")
    ax3[:set_ylabel]("Radius, r / m")

    for axes in fig[:get_axes]()
        axes[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 1))
    end

    fig[:tight_layout]()

    ax[:plot](x, T, linewidth=2)
    ax2[:plot](x, P, linewidth=2)
    ax3[:plot](x, r, linewidth=2)

    plt[:show]()
end


