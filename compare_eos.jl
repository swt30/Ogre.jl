using Ogre, PyCall, LaTeXStrings

import Ogre.plot
plot(a::Matrix, args...; kwargs...) = plot(a[:, 1], a[:, 2],
                                           args...; kwargs...)

function main()
    graph_h2o    = readcsv("data/seager-graphs/H2O.csv")
    graph_mgsio3 = readcsv("data/seager-graphs/MgSiO3.csv")
    graph_fe     = readcsv("data/seager-graphs/Fe.csv")

    pressures = logspace(7, 19)

    h2o          = map(Ogre.h2o_seager, pressures)
    mgsio3       = map(Ogre.mgsio3_seager, pressures)
    fe           = map(Ogre.fe_seager, pressures)
    my_h2o_300   = map(Ogre.my_h2o_300, pressures)
    my_h2o_500   = map(Ogre.my_h2o_500, pressures)
    my_h2o_800   = map(Ogre.my_h2o_800, pressures)
    my_h2o_1200  = map(Ogre.my_h2o_1200, pressures)

    plot(pressures, h2o, "k--", label="Seager H2O")
    plot(pressures, my_h2o_300, label="H2O (300K)")
    plot(pressures, my_h2o_500, label="H2O (500K)")
    plot(pressures, my_h2o_800, label="H2O (800K)")
    plot(pressures, my_h2o_1200, label="H2O (1200K)")

    ax = plt[:gca]()
    ax[:set_xlabel]("Pressure (P) / Pa")
    ax[:set_ylabel](L"Density ($\rho$) / kg m$^{-3}$")
    ax[:set_xlim]((1e7, 1e14))
    ax[:set_ylim]((1e2, 1e6))
    ax[:set_xscale]("log")
    ax[:set_yscale]("log")

    plt[:legend](loc=0)
    plt[:tight_layout]()
    plt[:show]()
end

main()
