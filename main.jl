using ogre, PyCall, LaTeXStrings
import Lazy: cycle
pygui()

# plot setup for Python
style = pyimport("matplotlib.style") # no @pyimport for now (Lint complains)
plt = pyimport("matplotlib.pyplot")
style[:use]("fivethirtyeight") # clean up this once dot-overloading is allowed

plot = plt[:plot]

# constants
const P_surface = 1.0e5
const mass_fractions = [1.]
const total_points = 100
const R_bracket = [0., 15.] .* R_earth

# this equation does not change with composition
const pressure_balance_eq = StructureEquation(pressure_balance)

# ...whereas the equation of mass continuity does
function setup_system(M::Real, eos::EOS)
    # planet layer options
    layer_densities = [eos]
    layer_edges = [0, cumsum(M.*mass_fractions)]

    # ODE system options
    m_inner, m_outer = layer_edges[1], layer_edges[end]
    solution_grid = linspace(m_outer, m_inner, total_points)

    # density-dependent equations change if layers or the EOS change
    density_profile = MassPiecewiseEOS(layer_densities, layer_edges)
    mass_continuity_with_eos(vs) = mass_continuity(vs, density_profile)
    mass_continuity_eq = StructureEquation(mass_continuity_with_eos)
    structure_equations = EquationSet([mass_continuity_eq,
                                       pressure_balance_eq])

    integrator.setup_find_radius(m_outer, mean(R_bracket), P_surface,
                                 structure_equations, solution_grid, R_bracket)
end

# radius finder for a solid sphere
function R(M::Real, eos::EOS; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1
    Rscale = in_earth_units ? 1/R_earth : 1

    planet_system = setup_system(M * Mscale, eos)
    r = find_radius!(planet_system) * Rscale

    r::Float64
end

# vectorized form of the above
function R{T<:Real}(ms::Vector{T}, eos::EOS; kwargs...)
    R_withEOS(M::T) = R(M::T, eos; kwargs...)
    map(R_withEOS, ms)
end

function main()
    ms = linspace(0.5, 10, 30) .* M_earth

    eoses = [eos.fe_seager,
             eos.h2o_seager,
             eos.mgsio3_seager]
    linestyles = cycle(["-"])
    colours = ["Crimson",
               "CornflowerBlue",
               "Sienna"]

    for (el, linestyle, colour) in zip(eoses, linestyles, colours)
        @time rs = R(ms, el)
        plt[:plot](ms ./ M_earth, rs ./ R_earth, label=el.fullname,
                   linestyle=linestyle, color=colour)
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
