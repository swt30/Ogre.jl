# STRUCTURE.JL
# Structural equations and planetary models

# CONSTANTS
#------------------------------------------------------------------------------

# Default values for the integrator
# TODO: re-document this module once permitted
module defaults
    import Ogre: R_earth

    # TODO: remove the hard coding on these
    "Surface pressure"
    const P_surf = 1.0e5
    "Surface temperature"
    const T_surf = 300
    "Mass fractions of the different layers"
    const mass_fractions = [1.]
    "Total number of points to use in the integration"
    const total_points = 100
    "Radius search bracket"
    const R_bracket = [0., 10.] .* R_earth
end

# STRUCTURAL EQUATIONS
#------------------------------------------------------------------------------

"A planetary structure equation that is not an EOS"
abstract StructureEquation{mc<:ModelComplexity} <: Equation

"Mass continuity: dr/dm = 1/4πr²ρ"
immutable MassContinuity{mc<:ModelComplexity, E<:EOS} <: StructureEquation{mc}
    eos::E
end
function MassContinuity{mc<:ModelComplexity}(eos::EOS{mc})
    MassContinuity{mc, typeof(eos)}(eos)
end

function Base.call(mce::MassContinuity, vs::ValueSet)
    dr_dm::Float64 = 0.0

    if isphysical(vs)
        r = radius(vs)
        ρ = mce.eos(vs)
        dr_dm = 1 / (4pi * r^2 * ρ)
    end

    dr_dm
end

"Pressure balance: dP/dm = -Gm/4πr⁴"
immutable PressureBalance <: StructureEquation{NoTemp} 
    # this equation does not change with composition
end
const pressurebalance = PressureBalance()

function Base.call(pbe::PressureBalance, vs::ValueSet)
    dP_dm::Float64 = 0.0

    if isphysical(vs)
        m = mass(vs)
        r = radius(vs)
        dP_dm = -(G * m) / (4pi * r^4)
    end

    dP_dm
end

"Adiabatic temperature gradient: dT/dm = -Gm/4πr⁴Cₚ"
immutable TemperatureGradient{E<:EOS, HC<:HeatCapacity} <: StructureEquation{WithTemp}
    eos::E
    heatcap::HC
end

function Base.call(tg::TemperatureGradient, pv::PhysicalValues)
    dT_dm::Float64 = 0.0

    if isphysical(pv)
        ρ = tg.eos(pv)
        cₚ = tg.heatcap(pv)
        r = radius(pv)
        m = mass(pv)
        T = temperature(pv)
        α = thermal_expansivity(pv, tg.eos)

        dT_dm = -(G * m * α * T) / (4pi * r^4 * ρ * cₚ)
    end

    dT_dm
end

using Calculus: derivative
partialx(F, x, y) = derivative(dx -> F(x + dx, y), 0.0)::Float64
partialy(F, x, y) = derivative(dy -> F(x, y + dy), 0.0)::Float64

"Thermal expansivity: αᵥ = -1/ρ (∂ρ/∂T)ₚ"
function thermal_expansivity(P, T, eos::EOS)
    FIX_THERMAL_EXPANSIVITY = false
    # TODO: remove this once exploratory data stuff is done

    if FIX_THERMAL_EXPANSIVITY
        α = 4e-4
    else
        dρdT = partialy(eos, P, T)
        ρ = eos(P, T)
        α = -1/ρ * dρdT
    end

    α::Float64
end

function thermal_expansivity(pv::PhysicalValues, eos::EOS)
    T = temperature(pv)
    P = pressure(pv)
    α = thermal_expansivity(P, T, eos)

    α::Float64
end

# PLANETARY STRUCTURE #
#------------------------------------------------------------------------------

typealias BoundaryValues ValueSet

"""Describes planetary parameters to be solved for an interior structure.

    * `M`: total mass
    * `structure_equations`: set of structural equations which incorporate the
      equation of state
    * `boundary_values`: a set of values specifying the external boundary
    * `solution_grid`: a grid of mass coordinates for the output
    * `radius_search_bracket`: the radius range to search when solving for R"""
abstract PlanetSystem{mc<:ModelComplexity}

"A planet with no temperature dependence"
immutable TempIndepPlanet{T<:Real, V<:AbstractVector} <: PlanetSystem{NoTemp}
    M::T
    structure_equations::EquationSet
    boundary_values::BoundaryValues{NoTemp}
    solution_grid::V
    radius_search_bracket::Vector{T}
end

"A planet with temperature dependence"
immutable TempDepPlanet{T<:Real, V<:AbstractVector} <: PlanetSystem{WithTemp}
    M::T
    structure_equations::EquationSet
    boundary_values::BoundaryValues{WithTemp}
    solution_grid::V
    radius_search_bracket::Vector{T}
end

function Base.show(io::IO, p::PlanetSystem)
    indepprefix = isa(p, TempIndepPlanet) ? "in" : ""
    dependent = indepprefix * "dependent"
    mass = p.M / M_earth
    g = p.solution_grid
    npoints = length(g)
    g1 = g[1] / M_earth
    gn = g[end] / M_earth
    rbracket = p.radius_search_bracket ./ R_earth
    
    println(io, "Temperature $dependent planet system with")
    println(io, "  mass: $mass M⊕")
    print(io, "  surface: ")
    show(io, p.boundary_values)
    println(io, "  $npoints mass points in $g1:$gn M⊕")
    println(io, "  radius bracket: $rbracket R⊕")
end

function PlanetSystem(M, eos::EOS{NoTemp}, bvs::BoundaryValues{NoTemp}, 
    grid=linspace(M, 0, defaults.total_points), r_bracket=defaults.R_bracket)

    masscontinuity = MassContinuity(eos)
    structure = EquationSet([masscontinuity, pressurebalance])

    TempIndepPlanet(M, structure, bvs, grid, r_bracket)
end
function PlanetSystem(M, eos::EOS{WithTemp}, Cₚ::HeatCapacity, 
    bvs::BoundaryValues{WithTemp}, grid=linspace(M, 0, defaults.total_points),
    r_bracket=defaults.R_bracket)

    masscontinuity = MassContinuity(eos)
    temperaturegradient = TemperatureGradient(eos, Cₚ)
    structure = EquationSet([masscontinuity, pressurebalance, temperaturegradient])

    TempDepPlanet(M, structure, bvs, grid, r_bracket)
end

"Default planet system (using values from module `defaults`)"
function DefaultPlanetSystem end
function DefaultPlanetSystem(M, eos::EOS{NoTemp})
    bvs = BoundaryValues(M, mean(defaults.R_bracket), defaults.P_surf)
    PlanetSystem(M, eos, bvs)
end
function DefaultPlanetSystem(M, eos::EOS{WithTemp}, Cₚ::HeatCapacity)
    bvs = BoundaryValues(M, mean(defaults.R_bracket), defaults.P_surf, defaults.T_surf)
    PlanetSystem(M, eos, Cₚ, bvs)
end

n_depvars{mc<:ModelComplexity}(sys::PlanetSystem{mc}) = n_depvars(mc)
npoints(sys::PlanetSystem) = length(sys.solution_grid)

"A planetary structure, containing mass grid `m` and internal physical values `y`"
abstract PlanetStructure{mc<:ModelComplexity}

"Planetary structure that holds only mass, radius, pressure"
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

"Planetary structure that holds temperature too"
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

"Generate a blank solution structure for a given planet system"
function blank_structure end
function blank_structure(sys::PlanetSystem{NoTemp})
    n = npoints(sys)
    MassRadiusPressureStructure(fill(NaN, 3, n))
end
function blank_structure(sys::PlanetSystem{WithTemp})
    n = npoints(sys)
    FullPlanetStructure(fill(NaN, 4, n))
end

"Current guess for planet radius, based on the search bracket"
R_guess(system::PlanetSystem) = mean(system.radius_search_bracket)

mass(ps::PlanetStructure) = sub(ps.data, 1, :)
radius(ps::PlanetStructure) = sub(ps.data, 2, :)
pressure(ps::PlanetStructure) = sub(ps.data, 3, :)
temperature(ps::PlanetStructure{WithTemp}) = sub(ps.data, 4, :)
gravity(ps::PlanetStructure) = G * mass(ps) ./ (radius(ps).^2)
nonmass(ps::PlanetStructure) = sub(ps.data, 2:nvars(ps), :)
function density(ps::PlanetStructure{NoTemp}, T, sys::PlanetSystem)
    eos = sys.structure_equations[1].eos
    map(P -> eos(P, T), pressure(ps))
end
function density(ps::PlanetStructure{WithTemp}, sys::PlanetSystem)
    eos = sys.structure_equations[1].eos
    map(eos, pressure(ps), temperature(ps))
end

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