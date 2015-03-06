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
abstract StructureEquation <: Equation

immutable MassContinuityEq <: StructureEquation
    equation::Function
end

immutable PressureBalanceEq <: StructureEquation
    equation::Function
end

immutable TemperatureGradientEq <: StructureEquation
    equation::Function
end

@doc "Mass continuity: dr/dm = 1/4πr²ρ" ->
function mass_continuity_f(vs::ValueSet, eos::EOS)
    dr_dm::Float64 = 0.0

    if isphysical(vs)
        r = radius(vs)
        ρ = eos(vs)
        dr_dm = 1 / (4pi * r^2 * ρ)
    end

    dr_dm
end
function MassContinuityEq(eos::EOS)
    mc(vs::ValueSet) = mass_continuity_f(vs, eos)
    MassContinuityEq(mc)
end


@doc "Pressure balance: dP/dm = -Gm/4πr⁴" ->
function pressure_balance_f(vs::ValueSet)
    dP_dm::Float64 = 0.0

    if isphysical(vs)
        m = mass(vs)
        r = radius(vs)
        dP_dm = -(G * m) / (4pi * r^4)
    end

    dP_dm
end
# this equation does not change with composition
PressureBalanceEq() = PressureBalanceEq(pressure_balance_f)
const pressurebalance = PressureBalanceEq()

@doc "Adiabatic energy gradient: dT/dm = -Gm/4πr⁴Cₚ" ->
function temperature_gradient_f(vs::PhysicalValues, eos::EOS,
    heatcap::HeatCapacity)

    dT_dm::Float64 = 0.0

    if isphysical(vs)
        ρ = eos(vs)
        Cₚ = heatcap(vs)
        r = radius(vs)
        m = mass(vs)

        dT_dm = -(G * m) / (4pi * r^4 * ρ * Cₚ)
    end

    dT_dm
end
function TemperatureGradientEq(eos::EOS, Cₚ::HeatCapacity)
    tg(vs::ValueSet) = temperature_gradient_f(vs, eos, Cₚ)
    TemperatureGradientEq(tg)
end


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
function PlanetSystem{T<:Real}(M::Real, eos::PressureEOS,
    bvs::MassRadiusPressure, grid::Vector{T}, r_bracket::Vector{T})

    masscontinuity = MassContinuityEq(eos)
    structure = EquationSet([masscontinuity, pressurebalance])
    PlanetSystem(M, structure, bvs, grid, r_bracket)
end
function PlanetSystem{T<:Real}(M::Real, eos::EOS, Cₚ::HeatCapacity,
    bvs::PhysicalValues, grid::Vector{T}, r_bracket::Vector{T})

    masscontinuity = MassContinuityEq(eos)
    temperaturegradient = TemperatureGradientEq(eos, Cₚ)
    structure = EquationSet([masscontinuity, pressurebalance, temperaturegradient])
    PlanetSystem(M, structure, bvs, grid, r_bracket)
end

@doc """
    A planetary structure, containing mass grid `m` and internal physical
    values `y`
    """ ->
type PlanetStructure{T<:Real}
    m::Vector{T}
    y::Matrix{T}
end

@doc "Generate a blank solution structure for a given planet system" ->
function blank_structure(sys::PlanetSystem)
    # array setup
    n_points = length(sys.solution_grid)
    t = fill(NaN, n_points)
    y = fill(NaN, (n_points, 2))
    solution = PlanetStructure(t, y)

    solution
end

@doc "Current guess for planet radius, based on the search bracket" ->
R_guess(system::PlanetSystem) = mean(system.radius_search_bracket)

# TODO: remove these funcs in favour of PlanetSystem constructor
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
    mass_continuity_with_eos(vs) = mass_continuity_f(vs, eos)
    masscontinuity = MassContinuityEq(mass_continuity_with_eos)
    structure_equations = EquationSet([masscontinuity,
                                       pressurebalance])

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


