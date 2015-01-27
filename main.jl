using Ogre, PyCall, LaTeXStrings
import Lazy: constantly

# plot setup for Python
pygui()
style = pyimport("matplotlib.style") # no @pyimport for now (Lint complains)
plt = pyimport("matplotlib.pyplot")
style[:use]("fivethirtyeight") # clean up this once dot-overloading is allowed
plot = plt[:plot]              # ditto

function main()
    ms = linspace(0.5, 10, 30) .* M_earth

    eoses = [Eos.fe_seager,
             Eos.h2o_seager,
             Eos.mgsio3_seager]
    linestyles = constantly("-")
    colours = ["Crimson", "CornflowerBlue", "Sienna"]

    for (el, ls, c) in zip(eoses, linestyles, colours)
        @time rs = R(ms, el)
        plot(ms ./ M_earth, rs ./ R_earth, label=el.fullname,
             linestyle=ls, color=c)
    end

    mh1 = readdlm("data/M-R/madhu/fe.out"; skipstart=1)
    mh2 = readdlm("data/M-R/madhu/perovskite.out"; skipstart=1)
    mh3 = readdlm("data/M-R/madhu/h2o.out"; skipstart=1)

    se1 = readcsv("data/M-R/seager/fe.csv")
    se2 = readcsv("data/M-R/seager/perovskite.csv")
    se3 = readcsv("data/M-R/seager/h2o.csv")

    plot(mh1[:, 1], mh1[:, 2], color="Black", linestyle=":")
    plot(mh2[:, 1], mh2[:, 2], color="Black", linestyle=":")
    plot(mh3[:, 1], mh3[:, 2], color="Black", linestyle=":",
         label="Madhu's curves")

    plot(se1[:, 1], se1[:, 2], color="Black", linestyle="--")
    plot(se2[:, 1], se2[:, 2], color="Black", linestyle="--")
    plot(se3[:, 1], se3[:, 2], color="Black", linestyle="--",
         label="Seager's curves")

    ax = plt[:gca]()
    ax[:set_xlabel](L"Mass / M$_\oplus$")
    ax[:set_ylabel](L"Radius / R$_\oplus$")
    ax[:set_xlim]((0, 10))
    ax[:set_ylim]((0, 3))
    ax[:set_xscale]("linear")
    ax[:set_yscale]("linear")
    xax = ax[:get_xaxis]()
    yax = ax[:get_yaxis]()

    plt[:legend](loc=0)

    ticker = pyimport("matplotlib.ticker")
    ScalarFormatter = ticker[:ScalarFormatter]
    xax[:set_major_formatter](ScalarFormatter())
    yax[:set_major_formatter](ScalarFormatter())
    plt[:tight_layout]()

    plt[:show]()
end

main()
