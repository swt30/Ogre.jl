module Integrator
using Ogre: Common, Constants, Structure, Eos
# Exported types
export BoundaryValues, PlanetSystem, PlanetStructure
# Exported functions
export zero, setup_system, find_radius, solve, R

#= Mass coordinates and other values =#

initial_values{T<:Real}(vs::ValueSet{T}) = [vs.r, vs.P]::Vector{T}
mass_coordinate(vs::ValueSet) = vs.m::Real

#= Types that define the planetary structure and the problem =#

typealias BoundaryValues ValueSet

immutable PlanetSystem{T<:Real}
    M::T
    structure_equations::EquationSet
    boundary_values::BoundaryValues
    solution_grid::Vector{T}
    radius_search_bracket::Vector{T}
end

type PlanetStructure{T<:Real}
    m::Vector{T}
    y::Matrix{T}
end

Base.zero(::Type{PlanetStructure}) = PlanetStructure([0.], [0. 0.])
R_guess(system::PlanetSystem) = mean(system.radius_search_bracket)

# this equation does not change with composition
const pressure_balance_eq = StructureEquation(pressure_balance)

# default values
const surface_pressure = 1.0e5
const mass_fractions = [1.]
const total_points = 100
const R_bracket = [0., 15.] .* R_earth

# set up a system with given EOS
# TODO: maybe shift this to be a constructor in the EOS module
function make_piecewise_EOS{T<:Real}(M::Real, eoses::Vector{SingleEOS},
    mass_fractions::Vector{T})

    layer_edges = [0, cumsum(M.*mass_fractions)]
    MassPiecewiseEOS(eoses, layer_edges) end

function setup_system{T<:SingleEOS, E<:Real}(M::Real, eoses::Vector{T},
    mass_fractions::Vector{E})

    eos = make_piecewise_EOS(M, eoses, mass_fractions)
    setup_system(M, eos)
end

function setup_system{T<:Real}(M::Real, eos::EOS,
    R_bracket::Vector{T}=R_bracket)

    # ODE system options
    m_inner, m_outer = 0, M
    solution_grid = linspace(m_outer, m_inner, total_points)

    # density-dependent equations change if layers or the EOS change
    mass_continuity_with_eos(vs) = mass_continuity(vs, eos)
    mass_continuity_eq = StructureEquation(mass_continuity_with_eos)
    structure_equations = EquationSet([mass_continuity_eq,
                                       pressure_balance_eq])

    setup_find_radius(m_outer, mean(R_bracket), surface_pressure,
                      structure_equations, solution_grid, R_bracket)
end

#= Functions to actually solve the structure, given boundary conditions =#

# generating blank solutions (filled with NaN)
function blank_structure(sys::PlanetSystem)
    # array setup
    n_points = length(sys.solution_grid)
    t = fill(NaN, n_points)
    y = fill(NaN, (n_points, 2))
    solution = PlanetStructure(t, y)

    solution
end

# solve and write to a given solution object
function solve!(sys::PlanetSystem, soln::PlanetStructure)
    # make an integrator and run it
    ode_func{T<:Real}(t::T, y::Vector{T}) = callfunc(sys.structure_equations,
                                                     t, y)
    x0 = initial_values(sys.boundary_values)
    t0 = mass_coordinate(sys.boundary_values)
    tgrid = sys.solution_grid

    solver = PlanetRK4(ode_func, x0, tgrid)
    for (i, x) in enumerate(solver)
        soln.m[i] = tgrid[i]
        soln.y[i, :] = x
    end
end

# generate a new solution object and solve
function solve(sys::PlanetSystem)
    soln = blank_structure(sys)
    solve!(sys, soln)

    soln
end

#= lower level functions to get a radius for a structure =#

# helper functions
notnan(value) = ~isnan(value)
dropnans(arr::Array) = filter(notnan, arr)
r_centre(result::PlanetStructure) = dropnans(result.y[:, 1])[end]
hit_the_centre(R::Real) = R < 0
not_far_enough(R::Real) = R > 100
hit_the_centre(result::PlanetStructure) = hit_the_centre(r_centre(result))
not_far_enough(result::PlanetStructure) = not_far_enough(r_centre(result))
unacceptable(result::PlanetStructure) = hit_the_centre(result) || not_far_enough(result)

function update_boundary_r(system::PlanetSystem)
    R_guess = mean(system.radius_search_bracket)
    updated_boundary_values = cpmod(system.boundary_values, r=R_guess)
    updated_system = cpmod(system, boundary_values=updated_boundary_values)
end

function update_R_search_bracket(system::PlanetSystem,
    r_low::Real, r_high::Real)

    cpmod(system, radius_search_bracket=[r_low, r_high])
end

function adapt_search_radius(system::PlanetSystem, result::PlanetStructure)
    if hit_the_centre(result)
        new_min_R = R_guess(system)
        new_max_R = system.radius_search_bracket[2]
        return update_R_search_bracket(system, new_min_R, new_max_R)
    elseif not_far_enough(result)
        new_min_R = system.radius_search_bracket[1]
        new_max_R = R_guess(system)
        return update_R_search_bracket(system, new_min_R, new_max_R)
    end
end

function find_radius(system::PlanetSystem)
    result = solve(system)
    system.radius_search_bracket
    if unacceptable(result)
        updated_search_radii = adapt_search_radius(system, result)
        prepped_for_re_solving = update_boundary_r(updated_search_radii)
        find_radius(prepped_for_re_solving)
    else
        system.boundary_values.r
    end
end

# set up a planetary system
function setup_find_radius{T<:Real}(M::T, R::T, P_surface::T,
    struct::EquationSet, solution_grid::Vector{T}, R_bracket::Vector{T})

    # boundary conditions and ODE setup
    bv = BoundaryValues(M, R, P_surface)
    system = PlanetSystem(M, struct, bv, solution_grid, R_bracket)
end

# create a system and find an appropriate radius
function find_radius{T<:Real}(M::T, structure::EquationSet, P_surface::T,
    solution_grid::Vector{T}, R_bracket::Vector{T})

    R_guess = mean(R_bracket)
    system = setup_system(M, R_guess, P_surface, structure,
                          solution_grid, R_bracket)
    radius = zero(PlanetStructure)
    R::T = find_radius(system)
end

#= higher level functions for doing MR diagrams =#

# radius finder for a solid sphere
function R(M::Real, eos::EOS; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1
    Rscale = in_earth_units ? 1/R_earth : 1

    planet_system = setup_system(M * Mscale, eos)
    r = find_radius(planet_system) * Rscale

    r::Float64
end

# vectorized form of the above
function R{T<:Real}(ms::Vector{T}, eos::EOS; in_earth_units=false)
    R_withEOS(M::T) = R(M::T, eos; in_earth_units=in_earth_units)
    map(R_withEOS, ms)
end

#= NUMERICAL METHODS =#

# accept either single or multi-valued starting conditions
typealias NumOrVec{T<:Real} Union(T, Vector{T})

# a single RK4 step: this function does the bulk of the work
function ode4_step{T<:Real}(F::Function, x::NumOrVec{T},
    tstart::T, tend::T)

    h::T = tend - tstart
    k = zeros(T, (length(x), 4))

    # Beginning derivative
    k[:,1] = h*F(tstart,        x)
    # Midstep derivatives
    k[:,2] = h*F(tstart + h./2, x + k[:,1]./2)
    k[:,3] = h*F(tstart + h./2, x + k[:,2]./2)
    # Ending derivative
    k[:,4] = h*F(tstart + h,    x + k[:,3])

    # Integrate
    x + k[:,1]./6 + k[:,2]./3 + k[:,3]./3 + k[:,4]./6
end

# types to handle solving
abstract IntegratorMethod
abstract FixedStepIntegrator <: IntegratorMethod
abstract RK4Integrator <: FixedStepIntegrator

typealias IntegratorState{I<:Integer, T<:Real} (I, NumOrVec{T})

# for general RK4 solving
immutable GenericRK4{T<:Real} <: RK4Integrator
    F::Function
    x0::NumOrVec{T}
    tgrid::Vector{T}
end

# for more specifically solving planetary structures
immutable PlanetRK4{T<:Real} <: RK4Integrator
    F::Function
    x0::NumOrVec{T}
    tgrid::Vector{T}
end

# integrator setup and steps
function Base.start(solver::FixedStepIntegrator)
    tindex = 1

    (tindex, solver.x0)
end

function Base.next(solver::RK4Integrator, state::IntegratorState)
    tindex, x = state

    # can't integrate past the end so just increment the t index instead
    if tindex == length(solver.tgrid)
        return x, (tindex + 1, x)
    end

    tstart = solver.tgrid[tindex]
    tend = solver.tgrid[tindex + 1]
    xn = ode4_step(solver.F, x, tstart, tend)

    newstate = (tindex + 1, xn)

    x, newstate
end

# termination conditions
function Base.done(solver::FixedStepIntegrator, state::IntegratorState)
    tindex, x = state
    tindex > length(solver.tgrid)
end

# in general, terminate when t<0 (m<0 for the case of planets). for planets, we
# let the integration continue *over* the centre and thus produce a zero value
# if r<0: this then signals the radius search function that we should increase
# our radius

# the length should be that of the solution grid
# (this allows us to use it in array comprehensions)
Base.length(solver::FixedStepIntegrator) = length(solver.tgrid)

# Wrappers to the above, providing dense output
function ode4_dense{T<:Real}(F::Function, x::T, tgrid::Vector{T})
    solver = GenericRK4(F, x, tgrid)

    # return a 1D array
    # the x[1] is to flatten any 1x1 arrays into single values
    [x[1] for x in solver]
end

function ode4_dense{T<:Real}(F::Function, x::Vector{T}, tgrid::Vector{T})
    solver = GenericRK4(F, x, tgrid)

    # return a 2D array
    reduce(hcat, solver)'
end

end
