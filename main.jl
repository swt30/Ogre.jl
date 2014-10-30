using ogre

# fixed parameters
const P_surface = 1e5
const layer_densities = [fe, mgsio3, h2o]
const mass_fractions = [1/3, 1/2, 1/6]
const total_points = 50
const R_bracket = [0., 5.] .* R_earth

function r(M)
    # planet layer options
    layer_edges = [0, cumsum(M.*mass_fractions)]

    # ODE system options
    m_inner, m_outer = layer_edges[1], layer_edges[end]
    solution_grid = linspace(m_outer, m_inner, total_points)

    # structure equations
    density_profile = PiecewiseEOS(layer_densities, layer_edges)
    mass_continuity_with_eos(vs) = mass_continuity(vs, density_profile)
    mass_continuity_eq = StructureEquation(mass_continuity_with_eos)
    pressure_balance_eq = StructureEquation(pressure_balance)
    structure_equations = EquationSet([mass_continuity_eq, pressure_balance_eq])

    # radius search
    get_radius(M, structure_equations, P_surface, solution_grid, R_bracket)
end

@vectorize_1arg Real r

ms = linspace(0.01, 10, 20) .* M_earth
@time rs = r(ms)

# python packages
using PyCall
pygui()
using PyPlot

# plot setup
@pyimport matplotlib.style as style
style.use("fivethirtyeight")

plot(ms ./ M_earth, rs ./ R_earth)
ax = gca()
ax[:set_xlabel](L"Mass / M$_\oplus$")
ax[:set_ylabel](L"Radius / R$_\oplus$")
ax[:set_xlim]((0, 10))
ax[:set_ylim]((0, 2))

tight_layout()
show()
