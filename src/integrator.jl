# INTEGRATOR.JL
# Numerical routines and solutions of planetary structral models

# Solutions
#------------------------------------------------------------------------------

@doc "Generate a blank solution structure for a given planet system" ->
function blank_structure(sys::PlanetSystem)
    # array setup
    n_points = length(sys.solution_grid)
    t = fill(NaN, n_points)
    y = fill(NaN, (n_points, 2))
    solution = PlanetStructure(t, y)

    solution
end

@doc """
    Solve a planetary structure, writing the solution to a given solution
    type.

    * `sys`: A `PlanetSystem` to solve
    * `soln`: A `PlanetStructure` to write the solution to
    """ ->
function solve!(sys::PlanetSystem, soln::PlanetStructure)
    # make an integrator and run it
    ode_func{T<:Real}(t::T, y::Vector{T}) = (sys.structure_equations(t, y))
    x0 = physical_values(sys.boundary_values)
    t0 = mass_coordinate(sys.boundary_values)
    tgrid = sys.solution_grid

    # put the values into the solution
    solver = PlanetRK4(ode_func, x0, tgrid)
    for (i, x) in enumerate(solver)
        soln.m[i] = tgrid[i]
        soln.y[i, :] = x
    end
end

@doc "Generate a new `PlanetStructure` object and solve a `PlanetSystem`" ->
function solve(sys::PlanetSystem)
    soln = blank_structure(sys)
    solve!(sys, soln)

    soln
end

# Radius searching
#------------------------------------------------------------------------------

# helper functions
notnan(value) = ~isnan(value)
dropnans(arr::Array) = filter(notnan, arr)
r_centre(result::PlanetStructure) = dropnans(result.y[:, 1])[end]
hit_the_centre(R::Real) = R < 0
not_far_enough(R::Real) = R > 100
hit_the_centre(result::PlanetStructure) = hit_the_centre(r_centre(result))
not_far_enough(result::PlanetStructure) = not_far_enough(r_centre(result))
unacceptable(result::PlanetStructure) = hit_the_centre(result) || not_far_enough(result)

@doc """
    Reproduce a `PlanetSystem` but change the radius value to be at the
    centre of its radius search bracket
    """ ->
function update_boundary_r(system::PlanetSystem)
    R_guess = mean(system.radius_search_bracket)
    updated_boundary_values = cpmod(system.boundary_values, r=R_guess)
    updated_system = cpmod(system, boundary_values=updated_boundary_values)
end

@doc "Reproduce a `PlanetSystem`, altering its radius search bracket" ->
function update_R_search_bracket(system::PlanetSystem,
    r_low::Real, r_high::Real)

    cpmod(system, radius_search_bracket=[r_low, r_high])
end

@doc """
    Based on a `PlanetStructure` solution for a given `PlanetSystem`,
    update the system with an adjusted radius search bracket
    """ ->
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

@doc """
    Recursively search the radius search bracket of a `PlanetSystem` to get a
    radius that produces an acceptable solution
    """ ->
function find_radius(system::PlanetSystem)
    result = solve(system)
    if unacceptable(result)
        updated_search_radii = adapt_search_radius(system, result)
        prepped_for_re_solving = update_boundary_r(updated_search_radii)
        find_radius(prepped_for_re_solving)
    else
        system.boundary_values.r
    end
end

@doc "Create and a solve for a `PlanetSystem`'s radius'" ->
function find_radius{T<:Real}(M::T, structure::EquationSet, P_surface::T,
    solution_grid::Vector{T}, R_bracket::Vector{T})

    R_guess = mean(R_bracket)
    system = setup_planet(M, R_guess, P_surface, structure,
                          solution_grid, R_bracket)
    radius = zero(PlanetStructure)
    R::T = find_radius(system)
end

# higher level functions for doing MR diagrams

@doc "Find the radius of a solid sphere of mass `M` using an given `EOS`." ->
function R(M::Real, eos::EOS; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1
    Rscale = in_earth_units ? 1/R_earth : 1

    planet_system = setup_planet(M * Mscale, eos)
    r = find_radius(planet_system) * Rscale

    r::Float64
end
# vectorized form of the above
function R{T<:Real}(ms::Vector{T}, eos::EOS; in_earth_units=false)
    R_withEOS(M::T) = R(M::T, eos; in_earth_units=in_earth_units)
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
    x + k[:,1]./6 + k[:,2]./3 + k[:,3]./3 + k[:,4]./6
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
# the integrator. We also define length so that we can put the results into
# array comprehensions, which need a length for the iterator in advance. This
# means that asking for the entire (dense) solution is really easy.

# iterator setup
function Base.start(solver::FixedStepIntegrator)
    tindex = 1

    (tindex, solver.x0)
end

# step function
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

# termination condition
function Base.done(solver::FixedStepIntegrator, state::IntegratorState)
    tindex, x = state
    # in general, terminate when we step off the end of the solution grid
    tindex > length(solver.tgrid)
end

# for planets, if r<0, the structure equations will yield zeroes: this allows
# us to step "over" the core, and so the negative sign of the radius at the
# central point signals the radius search function that we should increase our
# radius

# length will be that of the solution grid
Base.length(solver::FixedStepIntegrator) = length(solver.tgrid)

# Dense solvers
#------------------------------------------------------------------------------

@doc """
    Solve the ODE defined by function `F`, initial conds `x`, and fixed time
    steps `tgrid`. Returns dense output (an array of solutions at each point in
    `tgrid`)
""" ->
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
