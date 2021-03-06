# Structural equations and planetary models

using WaterData
import JLD, Dierckx, VoronoiDelaunay, GeometricalPredicates
import WaterData: istempdependent, extracteos


# General equation types

"Represents a physical equation"
abstract type Equation end
"Represents a set of equations which can be evaluated all at once"
struct EquationSet
    equations::Vector{Equation}
end

Base.length(es::EquationSet) = length(es.equations)
Base.getindex(es::EquationSet, i) = es.equations[i]
(es::EquationSet)(x::Real, y::Vector) = map(eq -> eq(x, y), es.equations)
(es::EquationSet)(vs::ValueSet) = map(eq -> eq(vs), es.equations)


# Structural equation types

"A planetary structure equation (that is not an EOS)"
abstract type StructureEquation <: Equation end

"Mass continuity: dr/dm = 1/4πr²ρ"
struct MassContinuity{E<:EOS} <: StructureEquation
    # this equation depends on density and therefore on the EOS
    eos::E
end

"Pressure balance: dP/dm = -Gm/4πr⁴"
struct PressureBalance <: StructureEquation
    # this equation does not change with composition
end
const pressurebalance = PressureBalance()

"Adiabatic temperature gradient: dT/dm = -GmαT/4πr⁴Cₚ"
struct TemperatureGradient <: StructureEquation
    # this equation depends on density and the heat capacity
    eos::EOS
    heatcap::HeatCapacity
end

"Thermal expansivity: αᵥ = -1/ρ (∂ρ/∂T)ₚ"
abstract type ThermalExpansivity <: StructureEquation end

struct GridThermalExp <: ThermalExpansivity
    alpha::EOS
end
struct ConstantThermalExp <: ThermalExpansivity
    alpha::Float64
end
struct NoThermalExp <: ThermalExpansivity; end

const water_thermexp = GridThermalExp(WaterData.load_full_eos()["thermexp"])
_thermexp(::EOS) = water_thermexp
_thermexp(::BME) = NoThermalExp()
_thermexp(::Vinet) = NoThermalExp()
function thermexp(eos::EOS, pv::PhysicalValues)
    singleeos = extracteos(eos, pv)
    _thermexp(singleeos)(pv)
end

# Evaluating structural equations

density(eos::EOS, vs) = eos(vs)

function (mce::MassContinuity)(vs::ValueSet)
    if isphysical(vs)
        r = radius(vs)
        P = pressure(vs)
        ρ = density(mce.eos, vs)
        return 1 / (4pi * r^2 * ρ)
    else
        return zero(Float64)
    end
end

function (::PressureBalance)(vs::ValueSet)
    if isphysical(vs)
        m = mass(vs)
        r = radius(vs)
        -(G * m) / (4pi * r^4)
    else
        zero(Float64)
    end
end

function (tg::TemperatureGradient)(pv::PhysicalValues)
    if isphysical(pv)
        ρ = tg.eos(pv)
        cₚ = tg.heatcap(pv)
        r = radius(pv)
        m = mass(pv)
        T = temperature(pv)
        if istempdependent(extracteos(tg.eos, pv))
            α = thermexp(tg.eos, pv)
            return -(G * m * α * T) / (4π * r^4 * ρ * cₚ)
        end
    end

    return zero(Float64)
end

(te::GridThermalExp)(pv::PhysicalValues) = te.alpha(pressure(pv), temperature(pv))
(te::ConstantThermalExp)(pv::PhysicalValues) = te.alpha(pressure(pv), temperature(pv))
(te::NoThermalExp)(pv::PhysicalValues) = te.alpha(pressure(pv), temperature(pv))

(te::GridThermalExp)(P::Pressure, T::Temperature) = te.alpha(P, T)
(te::ConstantThermalExp)(P::Pressure, T::Temperature) = te.alpha
(::NoThermalExp)(P::Pressure, T::Temperature) = zero(Float64)


# Planetary parameter types

const BoundaryValues = ValueSet

"""Describes planetary parameters to be solved for an interior structure.

* `M`: total mass
* `structure_equations`: set of structural equations which incorporate the
  equation of state
* `boundary_values`: a set of values specifying the external boundary
* `solution_grid`: a grid of mass coordinates for the output
* `radius_search_bracket`: the radius range to search when solving for R"""
abstract type PlanetSystem{mc<:ModelComplexity} end

"A planet with no temperature dependence"
mutable struct TempIndepPlanet <: PlanetSystem{NoTemp}
    M::Float64
    structure_equations::EquationSet
    boundary_values::BoundaryValues{NoTemp}
    solution_grid::Vector{Float64}
    radius_search_bracket::Vector{Float64}
end

"A planet with temperature dependence"
mutable struct TempDepPlanet <: PlanetSystem{WithTemp}
    M::Float64
    structure_equations::EquationSet
    boundary_values::BoundaryValues{WithTemp}
    solution_grid::Vector{Float64}
    radius_search_bracket::Vector{Float64}
    refine_surface!::Nullable{Function}
end

function PlanetSystem(M, eos::EOS, bvs::BoundaryValues{NoTemp},
                      grid=linspace(M, 0, defaults.total_points),
                      r_bracket=defaults.R_bracket)

    masscontinuity = MassContinuity(eos)
    structure = EquationSet([masscontinuity, pressurebalance])

    TempIndepPlanet(M, structure, bvs, grid, r_bracket)
end
function PlanetSystem(M, eos::EOS, Cₚ::HeatCapacity,
                      bvs::BoundaryValues{WithTemp},
                      grid=linspace(M, 0, defaults.total_points),
                      r_bracket=defaults.R_bracket, refine_surface=nothing)

    masscontinuity = MassContinuity(eos)
    temperaturegradient = TemperatureGradient(eos, Cₚ)
    structure = EquationSet([masscontinuity, pressurebalance, temperaturegradient])

    TempDepPlanet(M, structure, bvs, grid, r_bracket, refine_surface)
end

"Default planet system (using values from module `defaults`)"
function DefaultPlanetSystem end
function DefaultPlanetSystem(M, eos::EOS)
    bvs = BoundaryValues(M, mean(defaults.R_bracket), defaults.P_surf)
    PlanetSystem(M, eos, bvs)
end
function DefaultPlanetSystem(M, eos::EOS, Cₚ::HeatCapacity)
    bvs = BoundaryValues(M, mean(defaults.R_bracket), defaults.P_surf, defaults.T_surf)
    PlanetSystem(M, eos, Cₚ, bvs)
end


# Planet structural (solution) types

"A planetary structure, containing mass grid `m` and internal physical values `y`"
abstract type PlanetStructure{mc<:ModelComplexity} end

"Planetary structure that holds only mass, radius, pressure"
mutable struct MassRadiusPressureStructure <: PlanetStructure{NoTemp}
    data::Matrix{Float64}

    function MassRadiusPressureStructure(data::Matrix{Float64})
        if size(data)[1] == 3
            new(data)
        else
            error("Expected a size (3, n) array: got $(size(data))")
        end
    end
end
PlanetStructure(m, r, P) = MassRadiusPressureStructure(hcat(m, r, P)')

"Planetary structure that holds temperature too"
mutable struct FullPlanetStructure <: PlanetStructure{WithTemp}
    data::Matrix{Float64}

    function FullPlanetStructure(data::Matrix{Float64})
        if size(data)[1] == 4
            new(data)
        else
            error("Expected a size (4, n) array: got $(size(data))")
        end
    end
end
PlanetStructure(m, r, P, T) = FullPlanetStructure(hcat(m, r, P, T)')

"Generate a blank solution structure for a given planet system"
function blank_structure end
function blank_structure(sys::PlanetSystem{NoTemp})
    n = npoints(sys)
    m = nvars(NoTemp)
    MassRadiusPressureStructure(fill(NaN, m, n))
end
function blank_structure(sys::PlanetSystem{WithTemp})
    n = npoints(sys)
    m = nvars(WithTemp)
    FullPlanetStructure(fill(NaN, m, n))
end

# Interacting with planet parameter types

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

ndeps{mc<:ModelComplexity}(sys::PlanetSystem{mc}) = ndeps(mc)
npoints(sys::PlanetSystem) = length(sys.solution_grid)

maxradius(system) = system.radius_search_bracket[2]
minradius(system) = system.radius_search_bracket[1]
nextradiusguess(system) = mean(system.radius_search_bracket)
currentradiusguess(system) = system.boundary_values.r

function refine_r_bracket!(system::PlanetSystem, r_low, r_high)
    system.radius_search_bracket = [r_low, r_high]
end
"Choose a new radius guess for a planet system"
function refine_boundary_r!(system::PlanetSystem)
    system.boundary_values.r = nextradiusguess(system)
end
"Adjust the surface temperature/properties of a planet after each iteration"
function refine_surface!(system::TempDepPlanet)
    if !isnull(system.refine_surface!)
        do_refine! = get(system.refine_surface!)
        do_refine!(system.boundary_values)
    end
end
refine_surface!(system::TempIndepPlanet) = nothing

# Interacting with planet solution types

mass(ps::PlanetStructure) = ps.data[1, :]
setmass!(ps::PlanetStructure, i, m) = (ps.data[1, i] = m; nothing)
radius(ps::PlanetStructure) = ps.data[2, :]
setradius!(ps::PlanetStructure, i, r) = (ps.data[2, i] = r; nothing)
pressure(ps::PlanetStructure) = ps.data[3, :]
setpressure!(ps::PlanetStructure, i, P) = (ps.data[3, i] = P; nothing)
temperature(ps::PlanetStructure{WithTemp}) = ps.data[4, :]
settemperature!(ps::PlanetStructure, i, T) = (ps.data[4, i] = T; nothing)
gravity(ps::PlanetStructure) = G * mass(ps) ./ (radius(ps).^2)
nonmass(ps::PlanetStructure{NoTemp}) = ps.data[2:3, :]
nonmass(ps::PlanetStructure{WithTemp}) = ps.data[2:4, :]
setnonmass!(ps::PlanetStructure{NoTemp}, i, rP) = (ps.data[2:3, i] = rP; nothing)
setnonmass!(ps::PlanetStructure{WithTemp}, i, rPT) = (ps.data[2:4, i] = rPT; nothing)
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

centre(ps::PlanetStructure{NoTemp}) = MassRadiusPressure(ps.data[:, end]...)
centre(ps::PlanetStructure{WithTemp}) = PhysicalValues(ps.data[:, end]...)
surface(ps::PlanetStructure{NoTemp}) = MassRadiusPressure(ps.data[:, 1]...)
surface(ps::PlanetStructure{WithTemp}) = PhysicalValues(ps.data[:, 1]...)
