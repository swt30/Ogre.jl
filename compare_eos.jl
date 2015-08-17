using Ogre, PyPlot

type GraphValues
    P
    ρ
end
function GraphValues(filename::String)
    graph = readcsv(filename)
    P = graph[:, 1]
    ρ = graph[:, 2]
    GraphValues(P, ρ)
end

function main()
    graph_h2o    = GraphValues("data/eos/graphs/seager/h2o.csv")
    graph_mgsio3 = GraphValues("data/eos/graphs/seager/mgsio3.csv")
    graph_fe     = GraphValues("data/eos/graphs/seager/fe.csv")

    pressures    = logspace(5, 14)

    h2o          = map(Ogre.h2o_seager, pressures)
    mgsio3       = map(Ogre.mgsio3_seager, pressures)
    fe           = map(Ogre.fe_seager, pressures)
    my_h2o_300   = map(P -> Ogre.my_h2o_full(P, 300), pressures)
    my_h2o_500   = map(P -> Ogre.my_h2o_full(P, 500), pressures)
    my_h2o_800   = map(P -> Ogre.my_h2o_full(P, 800), pressures)
    my_h2o_1200  = map(P -> Ogre.my_h2o_full(P, 1200), pressures)

    fig = figure()
    ax = subplot(111)

    plot(pressures, h2o, "k--", label="BME/TFD (isothermal ice VII)")
    plot(pressures, my_h2o_300, label="Realistic water EOS (300K)")
    plot(pressures, my_h2o_500, label="(500K)")
    plot(pressures, my_h2o_800, label="(800K)")
    plot(pressures, my_h2o_1200, label="(1200K)")
    plot(graph_h2o.P, graph_h2o.ρ, label="Seager H2O")
    plot(graph_mgsio3.P, graph_mgsio3.ρ, label="Seager MgSiO3")
    plot(graph_fe.P, graph_fe.ρ, label="Seager Fe")

    xlabel("Pressure (P) / Pa")
    ylabel(L"Density ($\rho$) / kg m$^{-3}$")
    xlim(1e5, 1e14)
    # ylim(1e1, 1e5)
    xscale("log")
    yscale("log")

    legend(loc="best")
    tight_layout()
end

main()
