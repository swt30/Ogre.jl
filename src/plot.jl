module Plot
using Ogre: Common, Integrator

# plot setup for Python
using PyCall
pygui()
style = pyimport("matplotlib.style")
plt = pyimport("matplotlib.pyplot")
style[:use]("ggplot") # clean up this once dot-overloading is allowed

function plot(soln::Integrator.PlanetStructure)
    x = soln.m
    r = soln.y[:, 1]
    P = soln.y[:, 2]

    fig = plt[:figure]()
    ax = fig[:add_subplot](211)
    ax[:set_ylabel]("Pressure, P / Pa")

    ax2 = fig[:add_subplot](212, sharex=ax)
    ax2[:set_xlabel]("Mass enclosed, m / kg")
    ax2[:set_ylabel]("Radius, r / m")

    for axes in fig[:get_axes]()
        axes[:yaxis][:get_major_formatter]()[:set_powerlimits]((0, 1))
    end

    fig[:tight_layout]()

    ax[:plot](x, r, linewidth=2)
    ax2[:plot](x, P, linewidth=2)

    plt[:show]()
end


end