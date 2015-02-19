# STRUCTURE.JL
# Structural equations and planetary models

# CONSTANTS
#------------------------------------------------------------------------------

# default values for integrator
# TODO: remove the hard coding on these
const surface_pressure = 1.0e5
const mass_fractions = [1.]
const total_points = 100
const R_bracket = [0., 15.] .* R_earth

# STRUCTURAL EQUATIONS
#------------------------------------------------------------------------------

@doc "Type for a planetary structure equation that is not an EOS" ->
immutable StructureEquation <: Equation
    equation::Function
end

@doc "The mass continuity equation: dr/dm = 1/4πr²ρ" ->
function mass_continuity{T<:Real}(vs::ValueSet{T}, eos::EOS)
    if vs.m <= 0 || vs.r <= 0 || vs.P <= 0
        return 0.0
    end

    rho = eos(vs)
    dr_dm::Float64 = 1 / (4 * pi * vs.r^2 * rho)
end

@doc "The pressure balance equation: dP/dm = -Gm/4πr⁴" ->
function pressure_balance{T<:Real}(vs::ValueSet{T})
    if vs.m <= 0 || vs.r <= 0 || vs.P <= 0
        return 0.0
    end

    dP_dm::Float64 = -(G * vs.m) / (4 * pi * vs.r^4)
end

# this equation does not change with composition
const pressure_balance_eq = StructureEquation(pressure_balance)

# PLANETARY STRUCTURE #
#------------------------------------------------------------------------------

typealias BoundaryValues ValueSet

@doc """
    Describes planetary parameters to be solved for an interior structure.

    * `M`: total mass
    * `structure_equations`: set of structural equations which incorporate the
      equation of state
    * `boundary_values`: a set of values specifying the external boundary
    * `solution_grid`: a grid of mass coordinates for the output
    * `radius_search_bracket`: the radius range to search when solving for R
    """ ->
immutable PlanetSystem{T<:Real}
    M::T
    structure_equations::EquationSet
    boundary_values::BoundaryValues
    solution_grid::Vector{T}
    radius_search_bracket::Vector{T}
end

@doc """
    A planetary structure, containing mass grid `m` and internal physical
    values `y`
    """ ->
type PlanetStructure{T<:Real}
    m::Vector{T}
    y::Matrix{T}
end

Base.zero(::Type{PlanetStructure}) = PlanetStructure([0.], [0. 0.])

@doc "Current guess for planet radius, based on the search bracket" ->
R_guess(system::PlanetSystem) = mean(system.radius_search_bracket)

# set up a system with given EOS
@doc """
    Set up a `PlanetSystem` for radius finding.

    * `M`: total mass of the planet
    * `eos`: an equation of state (`EOS`) to be used for the mass continuity
      equation
    * `R_bracket` [optional]: Radius search bracket. Defaults to $R_bracket.
    """ ->
function setup_planet{T<:Real}(M::Real, eos::EOS,
    R_bracket::Vector{T}=R_bracket)

    # ODE system options
    m_inner, m_outer = 0, M
    solution_grid = linspace(m_outer, m_inner, total_points)

    # density-dependent equations change if layers or the EOS change
    mass_continuity_with_eos(vs) = mass_continuity(vs, eos)
    mass_continuity_eq = StructureEquation(mass_continuity_with_eos)
    structure_equations = EquationSet([mass_continuity_eq,
                                       pressure_balance_eq])

    setup_planet(m_outer, mean(R_bracket), surface_pressure,
                 structure_equations, solution_grid, R_bracket)
end
@doc """
    Lower-level function for setting up a planet when the structural
    equations have already been defined
    """ ->
function setup_planet{T<:Real}(M::T, R::T, P_surface::T, struct::EquationSet,
    solution_grid::Vector{T}, R_bracket::Vector{T})

    # boundary conditions and ODE setup
    bv = BoundaryValues(M, R, P_surface)
    system = PlanetSystem(M, struct, bv, solution_grid, R_bracket)
end
