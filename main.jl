using ogre, PyCall, PyPlot

# plot setup for Python
style = pyimport("matplotlib.style") # no @pyimport for now (Lint complains)
style[:use]("fivethirtyeight") # clean up this once dot-overloading is allowed

# constants
const P_surface = 1.0e5
const mass_fractions = [1.]
const total_points = 100
const R_bracket = [0., 15.] .* R_earth

# this equation does not change with composition
const pressure_balance_eq = StructureEquation(pressure_balance)

# fix some variables into the radius finding function
import ogre.integrator.get_radius
function get_radius(M::Real, structure_equations::EquationSet)
    get_radius(M, structure_equations, P_surface, solution_grid, R_bracket)
end

# radius finder for a solid sphere
function R(M::Real, eos::EOS)
    # planet layer options
    layer_densities = [eos]
    layer_edges = [0, cumsum(M.*mass_fractions)]

    # ODE system options
    m_inner, m_outer = layer_edges[1], layer_edges[end]
    solution_grid = linspace(m_outer, m_inner, total_points)

    # density-dependent equations change if layers or the EOS change
    density_profile = PiecewiseEOS(layer_densities, layer_edges)
    mass_continuity_with_eos(vs) = mass_continuity(vs, density_profile)
    mass_continuity_eq = StructureEquation(mass_continuity_with_eos)
    structure_equations = EquationSet([mass_continuity_eq,
                                       pressure_balance_eq])

    get_radius(M, structure_equations, P_surface, solution_grid, R_bracket)
end

# vectorized form of the above
function R{T<:Real}(ms::Vector{T}, eos::EOS)
    R_withEOS(M::T) = R(M::T, eos)
    @vectorize_1arg T R_withEOS
    R_withEOS(ms)
end

function main()
    ms = logspace(-1, 4, 50) .* M_earth

    for el in [h2o, mgsio3, fe]
        @time rs = R(ms, el)
        plot(ms ./ M_earth, rs ./ R_earth, label=el.name)
    end

    ax = gca()
    ax[:set_xlabel](L"Mass / M$_\oplus$")
    ax[:set_ylabel](L"Radius / R$_\oplus$")
    ax[:set_xlim]((0.1, 4000))
    ax[:set_ylim]((0.3, 11))
    ax[:set_xscale]("log")
    ax[:set_yscale]("log")
    xax = ax[:get_xaxis]()
    yax = ax[:get_yaxis]()

    legend(loc=0)

    ticker = pyimport("matplotlib.ticker")
    ScalarFormatter = ticker[:ScalarFormatter]
    xax[:set_major_formatter](ScalarFormatter())
    yax[:set_major_formatter](ScalarFormatter())
    tight_layout()

    show()
end

main()
