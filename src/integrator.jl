# INTEGRATOR.JL
# Numerical routines and solutions of planetary structral models

# Solutions
#------------------------------------------------------------------------------

@doc """Solve a `PlanetSystem`, writing the solution to a given
    `PlanetStructure`.

    * `sys`: A `PlanetSystem` to solve
    * `soln`: A `PlanetStructure` to write the solution to""" ->
function solve!(sys::PlanetSystem, soln::PlanetStructure)
    # make an integrator and run it
    ode_func{T<:Real}(t::T, y::Vector{T}) = (sys.structure_equations(t, y))
    x0 = nonmass(sys.boundary_values)
    t0 = mass(sys.boundary_values)
    tgrid = sys.solution_grid

    # put the values into the solution
    solver = PlanetRK4(ode_func, x0, tgrid)
    for (i, x) in enumerate(solver)
        mass(soln)[i] = tgrid[i]
        nonmass(soln)[:, i] = x
    end

    soln
end

@doc "Generate a new `PlanetStructure` object and solve a `PlanetSystem`" ->
function solve(sys::PlanetSystem)
    soln = blank_structure(sys)
    solve!(sys, soln)
end

# Radius searching
#------------------------------------------------------------------------------

# helper functions
notnan(value) = ~isnan(value)
hit_the_centre(R::Real) = R < 0
not_far_enough(R::Real) = R > 100
hit_the_centre(ps::PlanetStructure) = hit_the_centre(radius(centre(ps)))
not_far_enough(ps::PlanetStructure) = not_far_enough(radius(centre(ps)))
acceptable(ps::PlanetStructure) = !hit_the_centre(ps) && !not_far_enough(ps)

@doc """Reproduce a `PlanetSystem`, but change the radius guess to be at the
    centre of its radius search bracket""" ->
function update_boundary_r(system::PlanetSystem)
    R_guess = mean(system.radius_search_bracket)
    updated_boundary_values = cpmod(system.boundary_values, r=R_guess)
    updated_system = cpmod(system, boundary_values=updated_boundary_values)

    updated_system
end

@doc "Reproduce a `PlanetSystem`, altering its radius search bracket" ->
function update_R_search_bracket(system::PlanetSystem,
    r_low::Real, r_high::Real)

    cpmod(system, radius_search_bracket=[r_low, r_high])
end

@doc """Based on a `PlanetStructure` solution for a given `PlanetSystem`,
    update the system with an adjusted radius search bracket""" ->
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

@doc """Recursively search the radius search bracket of a `PlanetSystem` to
    get a result that has an acceptable error at the centre. Returns a
    `PlanetSystem` with the appropriate radius.""" ->
function converge(system::PlanetSystem)
    if R_guess(system) != system.boundary_values.r
        system = update_boundary_r(system)
    end
    result = solve(system)
    if acceptable(result)
        return system
    else
        system = adapt_search_radius(system, result)
        return converge(system)
    end
end

# basic functions for geting radii and structures
@doc "Get the correct internal structure for a `PlanetSystem`" ->
find_structure(system::PlanetSystem) = solve(converge(system))
@doc "Find the radius of a `PlanetSystem` by solving the system" ->
find_radius(system::PlanetSystem) = radius(surface(find_structure(system)))

@doc "Solve a `PlanetSystem` and give both the total radius and the structure" ->
function find_structure_and_radius(system::PlanetSystem)
    struct = find_structure(system)
    r = radius(surface(struct))

    return struct, r
end

# higher level functions for doing MR diagrams
@doc "Create and solve for a `PlanetSystem`'s radius'" ->
function find_radius{T<:Real}(M::T, structure::EquationSet, P_surface::T,
    solution_grid::Vector{T}, R_bracket::Vector{T})

    R_guess = mean(R_bracket)
    system = PlanetSystem(M, R_guess, P_surface, structure,
                          solution_grid, R_bracket)
    R::T = find_radius(system)
end

@doc "Find the radius of a solid sphere of mass `M` using an given `EOS`." ->
function R(M::Real, eos::EOS{NoTemp}; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1.
    Rscale = in_earth_units ? 1/R_earth : 1.

    sys = DefaultPlanetSystem(M * Mscale, eos)
    r = find_radius(sys) * Rscale

    r::Float64
end

function R(M::Real, eos::EOS{WithTemp}, Cₚ::HeatCapacity; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1.
    Rscale = in_earth_units ? 1/R_earth : 1.

    sys = DefaultPlanetSystem(M * Mscale, eos, Cₚ)
    r = find_radius(sys) * Rscale

    r::Float64
end

# vectorized form of the above
function R{T<:Real}(ms::Vector{T}, args...; in_earth_units=false)
    R_withEOS(M) = R(M, args...; in_earth_units=in_earth_units)
    map(R_withEOS, ms)
end

# Numerical methods
#------------------------------------------------------------------------------

# accept either single or multi-valued starting conditions
typealias NumOrVec{T<:Real} Union(T, Vector{T})

@doc "RK4 step of function `F`, initial conds `x`, from `tstart` to `tend`" ->
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
    # split onto multiple lines to help with code generation (?)
    (x + k[:,1]./6
       + k[:,2]./3
       + k[:,3]./3
       + k[:,4]./6)
end

# types to handle solving
abstract IntegratorMethod
abstract FixedStepIntegrator <: IntegratorMethod
abstract RK4Integrator <: FixedStepIntegrator

typealias IntegratorState{I<:Integer, T<:Real} (I, NumOrVec{T})

@doc "Type for general RK4 solving" ->
immutable GenericRK4{T<:Real} <: RK4Integrator
    F::Function
    x0::NumOrVec{T}
    tgrid::Vector{T}
end

@doc "Type for more specifically solving planetary structures" ->
immutable PlanetRK4{T<:Real} <: RK4Integrator
    F::Function
    x0::NumOrVec{T}
    tgrid::Vector{T}
end

# This section uses the iterator protocol by defining start, next, and done for
# the integrator. 

# iterator setup
function Base.start(solver::FixedStepIntegrator)
    (1, solver.x0)
end

# step function
function Base.next(solver::RK4Integrator, state::IntegratorState)
    tindex, x = state
    tnext = tindex + 1

    # can't integrate past the end so just increment the t index instead if we
    # hit the end of the solution grid
    if tindex == length(solver.tgrid)
        return x, (tnext, x)
    end

    # do a RK step to get the next values 
    tstart = solver.tgrid[tindex]
    tend = solver.tgrid[tnext]
    xnext = ode4_step(solver.F, x, tstart, tend)

    # the new state is the incremented index and the new values 
    newstate = (tnext, xnext)

    x, newstate
end

# termination condition
function Base.done(solver::FixedStepIntegrator, state::IntegratorState)
    tindex, x = state
    # we're done once we step off the end of the solution grid
    tindex > length(solver.tgrid) 
end

# for planets, if a value becomes unphysical (such as r<0 or P<0) then the
# structure equations will yield zeroes for each remaining step. This allows us
# to step "over" the core and leaves us with a slightly negative value at the
# centre. This negative sign signals the radius search function that we should
# increase our radius.

# We also define length so that we can put the results into array
# comprehensions, which need a length for the iterator in advance. This means
# that asking for the entire (dense) solution is really easy.
Base.length(solver::FixedStepIntegrator) = length(solver.tgrid)

# Dense solvers
#------------------------------------------------------------------------------

@doc """ Solve the ODE defined by function `F`, initial conds `x`, and fixed
    time steps `tgrid`. Returns dense output (an array of solutions at each
    point in `tgrid`) """ ->
function ode4_dense{T<:Real}(F::Function, x::T, tgrid::Vector{T})
    solver = GenericRK4(F, x, tgrid)

    # return a 1D array
    # the x[1] is to flatten any 1x1 arrays into single values
    [x[1] for x in solver]
end
function ode4_dense{T<:Real}(F::Function, x::Vector{T}, tgrid::Vector{T})
    solver = GenericRK4(F, x, tgrid)

    # return a 2D array
    hcat([x for x in solver]...)'
end
