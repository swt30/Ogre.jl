using Ogre, PyPlot

function mr_diagrams()
    ms = linspace(0.5, 10, 30) * M_earth

    Cₚ_const = Ogre.HeatCapacity(4200)  # constant heat capacity for liquid water
    Cₚ_inf = Ogre.HeatCapacity(Inf)  # infinite heat capacity
    Cₚ = Ogre.HeatCapacity("data/tabulated/heatcap-h2o.dat")  # realistic h.c.

    eoses = [Ogre.h2o_seager,
             Ogre.my_h2o_300,
             Ogre.my_h2o_500,
             Ogre.my_h2o_800]
    linestyles = ["--", "-", "-", "-"]
    colours = ["Black", "CornflowerBlue", "#fc4f30", "Gray"]
    legendtexts = ["BME/TFD (isothermal ice VII)",
                   "Realistic water EOS (isothermal 300K)",
                   "Realistic water EOS (isothermal 500K)",
                   "Realistic water EOS (isothermal 800K)"]

    fig = figure()
    ax = subplot(111)
    map(eoses, linestyles, colours, legendtexts) do el, ls, c, lg
        @time rs = Ogre.R(ms, el)
        plot(ms ./ M_earth, rs ./ R_earth, label=lg, linestyle=ls, color=c)
    end

    @time rs = Ogre.R(ms, Ogre.my_h2o_full, Cₚ)
    Ogre.plot(ms ./ M_earth, rs ./ R_earth, label="Realistic water EOS (adiabat)",
              linestyle="--", color="Crimson")

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

function phaseplots(Tsurf)
    M = M_earth
    R = R_earth

    eos = Ogre.my_h2o_full
    heatcap = Ogre.HeatCapacity("data/tabulated/heatcap-h2o.dat")
    Psurf = 1e5
    npoints = 200
    grid = linspace(M, 0, npoints)
    bcs = Ogre.ValueSet(M, R, Psurf, Tsurf)
    sys = Ogre.PlanetSystem(M, eos, heatcap, bcs, grid)

    sys = Ogre.converge(sys)
    soln = Ogre.solve(sys)
    close(:all)
    plot(soln, sys)
    tight_layout()
    transparent = (0,0,0,0)
    savefig("profiles-$(Tsurf)K.png", dpi=300, facecolor=transparent, edgecolor=nothing)
    Ogre.phaseplot(soln)
    savefig("phases-$(Tsurf)K.png", dpi=300, facecolor=transparent, edgecolor=nothing)
end

close(:all)
phaseplots(300)
# mr_diagrams()