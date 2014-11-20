using ogre, PyCall, PyPlot

# plot setup for Python
style = pyimport("matplotlib.style") # no @pyimport since Lint complains
style[:use]("ggplot") # clean up this once dot-overloading is allowed

function main()
    function mr_curve{T<:Real}(eos::EOS, ms::Vector{T})
        # fixed parameters
        const P_surface = 1.0e5
        const layer_densities = [eos]
        const mass_fractions = [1.]
        const total_points = 100
        const R_bracket = [0., 15.] .* R_earth

        function r(M::T)
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
            structure_equations = EquationSet([mass_continuity_eq,
                                               pressure_balance_eq])

            # radius search
            get_radius(M, structure_equations, P_surface, solution_grid, R_bracket)
        end

        @vectorize_1arg Real r

        rs = r(ms)
    end

    ms = logspace(-1, 4, 50) .* M_earth

    for el in [h2o, mgsio3, fe]
        @time rs = mr_curve(el, ms)
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

    legend()

    ticker = pyimport("matplotlib.ticker")
    ScalarFormatter = ticker[:ScalarFormatter]
    xax[:set_major_formatter](ScalarFormatter())
    yax[:set_major_formatter](ScalarFormatter())
    tight_layout()

    show()
end

main()
