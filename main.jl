using Ogre, LaTeXStrings
import Iterators: repeated

using Dierckx
using PyPlot

function main()
    ms = linspace(0.5, 10, 30) * M_earth

    Cₚ_func(T::Real) = 4200
    Cₚ = Ogre.HeatCapacity(Cₚ_func)
    Cₚ_inf_func(T::Real) = Inf
    Cₚ_inf = Ogre.HeatCapacity(Cₚ_inf_func)

    eoses = [Ogre.h2o_seager,
             Ogre.my_h2o_300,
             Ogre.my_h2o_500]
    linestyles = ["--", "-", "-"]
    colours = ["Black", "CornflowerBlue", "#fc4f30"]

    Ogre.R(M_earth, eoses[1])

    legendtexts = ["BME/TFD (isothermal ice VII)",
                   "Realistic water EOS (isothermal 300K)",
                   "Realistic water EOS (isothermal 500K)"]

    fig = figure()
    map(eoses, linestyles, colours, legendtexts) do el, ls, c, lg
        @time rs = Ogre.R(ms, el)
        plot(ms ./ M_earth, rs ./ R_earth, label=lg,
             linestyle=ls, color=c)
    end
    # @time rs = Ogre.R(ms, Ogre.my_h2o, Cₚ)
    # Ogre.plot(ms ./ M_earth, rs ./ R_earth, label="Realistic water EOS (adiabat)",
    #      linestyle="-", color="Crimson")

    ax = gca()
    xlabel(L"Mass / M$_\oplus$")
    ylabel(L"Radius / R$_\oplus$")
    xlim(0, 10)
    ylim(0, 5)
    xscale("linear")
    yscale("linear")
    xax = ax[:get_xaxis]()
    yax = ax[:get_yaxis]()

    legend(loc=0)

    ScalarFormatter = ticker.ScalarFormatter
    xax[:set_major_formatter](ScalarFormatter())
    yax[:set_major_formatter](ScalarFormatter())

    tight_layout()

end
