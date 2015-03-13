using Ogre, LaTeXStrings
import Iterators: repeated

using Dierckx

function main()
    ms = linspace(0.5, 10, 30) * M_earth

    Cₚ_func(T::Real) = 4200
    Cₚ = Ogre.HeatCapacity(Cₚ_func)
    Cₚ_inf_func(T::Real) = Inf
    Cₚ_inf = Ogre.HeatCapacity(Cₚ_inf_func)

    eoses = [Ogre.h2o_seager,
             Ogre.my_h2o_300]
    linestyles = repeated("-")
    colours = ["Crimson", "CornflowerBlue", "Sienna",
                "LimeGreen", "LightSteelBlue"]

    Ogre.R(M_earth, eoses[1])

    map(eoses, linestyles, colours) do el, ls, c
        @time rs = Ogre.R(ms, el)
        Ogre.plot(ms ./ M_earth, rs ./ R_earth, label=el.fullname,
             linestyle=ls, color=c)
    end
    @time rs = Ogre.R(ms, Ogre.my_h2o, Cₚ)
    Ogre.plot(ms ./ M_earth, rs ./ R_earth, label=Ogre.my_h2o.fullname,
         linestyle="--", color="Black")

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

    ScalarFormatter = ticker[:ScalarFormatter]
    xax[:set_major_formatter](ScalarFormatter())
    yax[:set_major_formatter](ScalarFormatter())
    plt[:tight_layout]()

    plt[:show]()
end
