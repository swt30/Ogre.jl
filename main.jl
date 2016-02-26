using Ogre
using WaterData
using PyPlot
using Colors
using DataStructures: OrderedDict

function mr_diagrams()
    ms = linspace(0.5M_earth, 10M_earth, 3)

    h2o = WaterData.load_full_eos()["grid"]
    funcs = WaterData.load_piecewise_eoses()
    mgsio3 = funcs["mgsio3"]
    fe = funcs["fe"]

    eoses = [fe, mgsio3, h2o]
    mass_fractions = [0.15, 0.30, 0.55]
    total_points = 200
    R_bracket = [0, 10] * R_earth

    # Cₚ_const = WaterData.ConstantHeatCapacity(4200)  # constant heat capacity for liquid water
    # Cₚ_inf = WaterData.ConstantHeatCapacity(Inf)  # infinite heat capacity
    Cₚ = WaterData.load_heat_capacity()["heatcap_h2o"] # realistic h.c.

    isotherms = OrderedDict(# 300 => WaterData.slice(h2o, 300),
                            400 => WaterData.slice(h2o, 400),
                            500 => WaterData.slice(h2o, 500),
                            600 => WaterData.slice(h2o, 600),
                            700 => WaterData.slice(h2o, 700),
                            800 => WaterData.slice(h2o, 800),
                            # 900 => WaterData.slice(h2o, 900),
                            # 1000 => WaterData.slice(h2o, 1000)
                            )

    function R_adiabat(M, Tsurf, Psurf)
        Me = M / M_earth
        bvs = Ogre.ValueSet(M, R_earth, Psurf, Tsurf)
        grid = linspace(M, 0, total_points)
        planet = Ogre.PlanetSystem(M, h2o, Cₚ, bvs, grid, R_bracket)

        Ogre.find_radius!(planet)
    end

    function R_adiabat_cored(M, Tsurf, Psurf)
        Me = M / M_earth
        bvs = Ogre.ValueSet(M, R_earth, Psurf, Tsurf)
        grid = linspace(M, 0, total_points)
        coredeos = Ogre.MassPiecewiseEOS(eoses, M, mass_fractions)
        planet = Ogre.PlanetSystem(M, coredeos, Cₚ, bvs, grid, R_bracket)

        Ogre.find_radius!(planet)
    end

    function R_isotherm(M, Tsurf, Psurf)
        bvs = Ogre.ValueSet(M, R_earth, Psurf)
        grid = linspace(M, 0, total_points)
        planet = Ogre.PlanetSystem(M, isotherms[Tsurf], bvs, grid, R_bracket)

        Ogre.find_radius!(planet)
    end

    temps = keys(isotherms)
    dropfirst(x) = drop(x, 1)

    for f in [R_isotherm, R_adiabat, R_adiabat_cored], Psurf in [1e5, 1e6, 1e7]
        cmap = Colors.colormap("Reds", length(temps) + 1) |> dropfirst
        figure()
        for (T, color) in zip(temps, cmap)
            c = map(f -> f(color), [red, green, blue])
            @time rs = map(M -> f(M, T, Psurf), ms)
            plot(ms ./ M_earth, rs ./ R_earth, color=c, label="$T K")
        end

        method = string(f)[3:end]
        log10P = floor(Int, log10(Psurf))
        title("Psurf=1e$log10P Pa; $method")
        xlabel(L"Mass / M$_⊕$")
        ylabel(L"Radius / R$_⊕$")
        legend(loc="best")
        xlim((0, 10))
        ylim((1, 5))
        tight_layout()
        savefig("$(method)_$(log10P)Pa.png")
        close()
    end
end

function phaseplots(M, Psurf, Tsurf)
    M = M * M_earth
    R = cbrt(M)
    npoints = 200
    grid = linspace(M, 0, npoints)
    eos = WaterData.load_full_eos()["grid"]
    heatcap = WaterData.load_heat_capacity()["heatcap_h2o"]

    bcs = Ogre.ValueSet(M, R, Psurf, Tsurf)
    sys = Ogre.PlanetSystem(M, eos, heatcap, bcs, grid)

    soln = Ogre.find_structure!(sys)
    plot(soln, sys)
    # tight_layout()
    Ogre.phaseplot(soln)
end

close(:all)
phaseplots(1, 1e6, 400)
# mr_diagrams()
