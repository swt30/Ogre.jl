# STRUCTURE.JL
# Structural equations and planetary models

# CONSTANTS
#------------------------------------------------------------------------------

# default values for integrator
# TODO: remove the hard coding on these
const P_surf = 1.0e5
const T_surf = 300
const mass_fractions = [1.]
const total_points = 100
const R_bracket = [0., 10.] .* R_earth

# STRUCTURAL EQUATIONS
#------------------------------------------------------------------------------

@doc "Type for a planetary structure equation that is not an EOS" ->
abstract StructureEquation{mc<:ModelComplexity} <: Equation

immutable MassContinuityEq{mc<:ModelComplexity} <: StructureEquation{mc}
    equation::Function
end

immutable PressureBalanceEq <: StructureEquation{NoTemp}
    equation::Function
end

immutable TemperatureGradientEq <: StructureEquation{WithTemp}
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
function MassContinuityEq{mc<:ModelComplexity}(eos::EOS{mc})
    masscontinuity(vs::ValueSet) = mass_continuity_f(vs, eos)
    MassContinuityEq{mc}(masscontinuity)
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
    tempgradient(vs::ValueSet) = temperature_gradient_f(vs, eos, Cₚ)
    TemperatureGradientEq(tempgradient)
end

# PLANETARY STRUCTURE #
#------------------------------------------------------------------------------

typealias BoundaryValues ValueSet

@doc """Describes planetary parameters to be solved for an interior structure.

    * `M`: total mass
    * `structure_equations`: set of structural equations which incorporate the
      equation of state
    * `boundary_values`: a set of values specifying the external boundary
    * `solution_grid`: a grid of mass coordinates for the output
    * `radius_search_bracket`: the radius range to search when solving for R
    """ ->
abstract PlanetSystem{mc<:ModelComplexity}

immutable TempIndepPlanet{T<:Real} <: PlanetSystem{NoTemp}
    M::T
    structure_equations::EquationSet
    boundary_values::BoundaryValues{NoTemp}
    solution_grid::Vector{T}
    radius_search_bracket::Vector{T}
end
immutable TempDepPlanet{T<:Real} <: PlanetSystem{WithTemp}
    M::T
    structure_equations::EquationSet
    boundary_values::BoundaryValues{WithTemp}
    solution_grid::Vector{T}
    radius_search_bracket::Vector{T}
end

function PlanetSystem{T<:Real}(M::T, eos::EOS{NoTemp},
    bvs::BoundaryValues{NoTemp}, grid::Vector{T}=linspace(M, 0, total_points),
    r_bracket::Vector{T}=R_bracket)

    masscontinuity = MassContinuityEq(eos)
    structure = EquationSet([masscontinuity, pressurebalance])

    TempIndepPlanet(M, structure, bvs, grid, r_bracket)
end
function PlanetSystem{T<:Real}(M::T, eos::EOS{WithTemp}, Cₚ::HeatCapacity,
    bvs::BoundaryValues{WithTemp}, grid::Vector{T}=linspace(M, 0, total_points),
    r_bracket::Vector{T}=R_bracket)

    masscontinuity = MassContinuityEq(eos)
    temperaturegradient = TemperatureGradientEq(eos, Cₚ)
    structure = EquationSet([masscontinuity, pressurebalance, temperaturegradient])

    TempDepPlanet(M, structure, bvs, grid, r_bracket)
end

function DefaultPlanetSystem(M::Real, eos::EOS{NoTemp})
    bvs = BoundaryValues(M, mean(R_bracket), P_surf)
    PlanetSystem(M, eos, bvs)
end
function DefaultPlanetSystem(M::Real, eos::EOS{WithTemp}, Cₚ::HeatCapacity)
    bvs = BoundaryValues(M, mean(R_bracket), P_surf, T_surf)
    PlanetSystem(M, eos, Cₚ, bvs)
end

n_depvars{mc<:ModelComplexity}(sys::PlanetSystem{mc}) = n_depvars(mc)
npoints(sys::PlanetSystem) = length(sys.solution_grid)

@doc """A planetary structure, containing mass grid `m` and internal physical
    values `y` """ ->
abstract PlanetStructure{mc<:ModelComplexity}

type MassRadiusPressureStructure{T<:Real} <: PlanetStructure{NoTemp}
    data::Matrix{T}

    function MassRadiusPressureStructure(data::Matrix{T})
        if size(data)[1] == 3
            new(data)
        else
            error("Expected a size (3, n) array: got $(size(data))")
        end
    end
end
function MassRadiusPressureStructure{T<:Real}(data::Matrix{T})
    MassRadiusPressureStructure{T}(data)
end
PlanetStructure(m, r, P) = MassRadiusPressureStructure(hcat(m, r, P)')

type FullPlanetStructure{T<:Real} <: PlanetStructure{WithTemp}
    data::Matrix{T}

    function FullPlanetStructure(data::Matrix{T})
        if size(data)[1] == 4
            new(data)
        else
            error("Expected a size (4, n) array: got $(size(data))")
        end
    end
end
function FullPlanetStructure{T<:Real}(data::Matrix{T})
    FullPlanetStructure{T}(data)
end
PlanetStructure(m, r, P, T) = FullPlanetStructure(hcat(m, r, P, T)')

@doc "Generate a blank solution structure for a given planet system" ->
function blank_structure(sys::PlanetSystem{NoTemp})
    n = npoints(sys)
    MassRadiusPressureStructure(fill(NaN, 3, n))
end
function blank_structure(sys::PlanetSystem{WithTemp})
    n = npoints(sys)
    FullPlanetStructure(fill(NaN, 4, n))
end

@doc "Current guess for planet radius, based on the search bracket" ->
R_guess(system::PlanetSystem) = mean(system.radius_search_bracket)

mass(ps::PlanetStructure) = sub(ps.data, 1, :)
radius(ps::PlanetStructure) = sub(ps.data, 2, :)
pressure(ps::PlanetStructure) = sub(ps.data, 3, :)
temperature(ps::PlanetStructure{WithTemp}) = sub(ps.data, 4, :)

npoints(ps::PlanetStructure) = length(mass(ps))
ndeps{MC<:ModelComplexity}(::PlanetStructure{MC}) = ndeps(MC)
nvars{MC<:ModelComplexity}(::PlanetStructure{MC}) = nvars(MC)

function centre(ps::PlanetStructure)
    centrecoords = (:, npoints(ps))
    ValueSet(sub(ps.data, centrecoords)...)
end

function surface(ps::PlanetStructure)
    surfcoords = (:, 1)
    ValueSet(sub(ps.data, surfcoords)...)
end

# some convenience functions
nonmass(ps::PlanetStructure) = sub(ps.data, 2:nvars(ps), :)