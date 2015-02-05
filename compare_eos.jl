using Ogre: Eos
using PyCall, LaTeXStrings

style = pyimport("matplotlib.style")
plt = pyimport("matplotlib.pyplot")
style[:use]("fivethirtyeight")

plot(a::Matrix, args...; kwargs...) = plot(a[:, 1], a[:, 2],
                                           args...; kwargs...)
plot(args...; kwargs...) = plt[:plot](args...; kwargs...)


function main()
    graph_h2o    = readcsv("data/seager-graphs/H2O.csv")
    graph_mgsio3 = readcsv("data/seager-graphs/MgSiO3.csv")
    graph_fe     = readcsv("data/seager-graphs/Fe.csv")

    pressures = logspace(7, 19)

    h2o          = map(P -> callfunc(Eos.h2o_seager, P), pressures)
    mgsio3       = map(P -> callfunc(Eos.mgsio3_seager, P), pressures)
    fe           = map(P -> callfunc(Eos.fe_seager, P), pressures)

    plot(pressures, h2o, label="H2O")
    plot(pressures, mgsio3, label="MgSiO3")
    plot(pressures, fe, label="Fe")
    plot(graph_h2o, "k-", linewidth=1, label="Seager's functions")
    plot(graph_mgsio3, "k-", linewidth=1)
    plot(graph_fe, "k-", linewidth=1)

    ax = plt[:gca]()
    ax[:set_xlabel]("Pressure (P) / Pa")
    ax[:set_ylabel](L"Density ($\rho$) / kg m$^{-3}$")
    ax[:set_xlim]((1e7, 1e19))
    ax[:set_ylim]((1e2, 1e8))
    ax[:set_xscale]("log")
    ax[:set_yscale]("log")

    plt[:legend](loc=0)
    plt[:tight_layout]()
    plt[:show]()
end

main()
