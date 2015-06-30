# INTEGRATOR.JL
# Numerical routines for solutions of planetary structure models

using Compat

# Solutions
#------------------------------------------------------------------------------

"Solve a `PlanetSystem`, writing the solution to a given `PlanetStructure`."
function solve!(sys::PlanetSystem, soln::PlanetStructure)
    # initial values
    x0 = nonmass(sys.boundary_values)
    t0 = mass(sys.boundary_values)
    tgrid = sys.solution_grid

    # solve and write to solution array
    solver = PlanetRK4(sys.structure_equations, x0, tgrid)
    for (i, x) in enumerate(solver)
        # mass and nonmass are array views, so we can assign to them 
        mass(soln)[i] = tgrid[i]
        nonmass(soln)[:, i] = x
    end

    soln
end

"Generate a new `PlanetStructure` object and solve a `PlanetSystem`"
function solve(sys::PlanetSystem)
    soln = blank_structure(sys)
    solve!(sys, soln)
end

# Radius searching
#------------------------------------------------------------------------------

# helper functions
"Has the radius passed the central point?"
function hit_the_centre end
hit_the_centre(R::Real) = R < 0
hit_the_centre(ps::PlanetStructure) = hit_the_centre(radius(centre(ps)))
"Is the radius not yet close enough to the centre?"
function not_far_enough end
not_far_enough(R::Real) = R > 100
not_far_enough(ps::PlanetStructure) = not_far_enough(radius(centre(ps)))
"Is the structural solution acceptable?"
acceptable(ps::PlanetStructure) = !hit_the_centre(ps) && !not_far_enough(ps)

""" Clone a `PlanetSystem`, updating the radius guess appropriately

    The new system will have R_guess at the centre of its radius search
    bracket. """
function update_boundary_r(system::PlanetSystem)
    R_guess = mean(system.radius_search_bracket)
    updated_boundary_values = cpmod(system.boundary_values, r=R_guess)
    updated_system = cpmod(system, boundary_values=updated_boundary_values)

    updated_system
end

"Clone a `PlanetSystem`, updating its radius search bracket"
function update_R_search_bracket(system::PlanetSystem,
    r_low::Real, r_high::Real)

    cpmod(system, radius_search_bracket=[r_low, r_high])
end

""" Update a `PlanetSystem` with an appropriate search radius 

    The radius is determined based on the previous `PlanetStructure` solution
    """
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

""" Repeatedly solve a `PlanetSystem` until its radius is suitable

    Recursively search the radius search bracket of a `PlanetSystem` to
    get a result that has an acceptable error at the centre. Returns a
    `PlanetSystem` with the appropriate radius. """
function converge(system::PlanetSystem)
    if R_guess(system) ≠ system.boundary_values.r
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
"Get the correct internal structure for a `PlanetSystem`"
find_structure(system::PlanetSystem) = solve(converge(system))
"Find the radius of a `PlanetSystem` by solving the system"
find_radius(system::PlanetSystem) = radius(surface(find_structure(system)))

"Solve a `PlanetSystem` and give both the total radius and the structure"
function find_structure_and_radius(system::PlanetSystem)
    struct = find_structure(system)
    r = radius(surface(struct))

    return struct, r
end

# higher level functions for doing MR diagrams
"Create and solve for a `PlanetSystem`'s radius'"
function find_radius{T<:Real}(M::T, structure::EquationSet, P_surface::T,
    solution_grid::Vector{T}, R_bracket::Vector{T})

    R_guess = mean(R_bracket)
    system = PlanetSystem(M, R_guess, P_surface, structure,
                          solution_grid, R_bracket)
    R::T = find_radius(system)
end

"Find the radius of a solid sphere of mass `M` using an given `EOS`."
function R(M::Real, eos::EOS{NoTemp}; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1.
    Rscale = in_earth_units ? 1/R_earth : 1.

    sys = DefaultPlanetSystem(M * Mscale, eos)
    r = find_radius(sys) * Rscale

    r::Float64
end

# temperature-dependent version
function R(M::Real, eos::EOS{WithTemp}, Cₚ::HeatCapacity; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1.
    Rscale = in_earth_units ? 1/R_earth : 1.

    sys = DefaultPlanetSystem(M * Mscale, eos, Cₚ)
    r = find_radius(sys) * Rscale

    r::Float64
end

# vectorized form of the above
function R{T<:Real}(ms::AbstractVector{T}, args...; in_earth_units=false)
    R_withEOS(M) = R(M, args...; in_earth_units=in_earth_units)
    map(R_withEOS, ms)
end

# Numerical methods
#------------------------------------------------------------------------------

# accept either single or multi-valued starting conditions
typealias NumOrVec{T<:Real} Union(T, Vector{T})

"Do a RK4 step of function `F`, initial conds `x`, from `tstart` to `tend`"
function ode4_step(F, x, tstart, tend)
    h = tend - tstart
    k = zeros(eltype(x), (length(x), 4))

    # Beginning derivative
    k[:,1] = h*F(tstart,        x)
    # Midstep derivatives
    k[:,2] = h*F(tstart + h./2, x + k[:,1]./2)
    k[:,3] = h*F(tstart + h./2, x + k[:,2]./2)
    # Ending derivative
    k[:,4] = h*F(tstart + h,    x + k[:,3])

    # Integrate
    (x + k[:,1]/6 
       + k[:,2]/3 
       + k[:,3]/3 
       + k[:,4]/6)
end

# types to handle solving
abstract IntegratorMethod
abstract FixedStepIntegrator <: IntegratorMethod
abstract RK4Integrator <: FixedStepIntegrator

"General RK4 solver"
immutable GenericRK4{T<:Real, V<:AbstractVector{Float64}} <: RK4Integrator
    func
    x0::NumOrVec{T}
    tgrid::V
end

"Planetary structure RK4 solver"
immutable PlanetRK4{T<:Real, V<:AbstractVector{Float64}} <: RK4Integrator
    func
    x0::NumOrVec{T}
    tgrid::V
end

# This section uses the iterator protocol by defining start, next, and done for
# the integrator. 

# iterator setup
function Base.start(solver::FixedStepIntegrator)
    (1, solver.x0)
end

# step function
function Base.next(solver::RK4Integrator, state)
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
    xnext = ode4_step(solver.func, x, tstart, tend)

    # the new state is the incremented index and the new values 
    newstate = (tnext, xnext)

    x, newstate
end

# termination condition
function Base.done(solver::FixedStepIntegrator, state)
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

""" Solve the ODE dx/dt = `F` with initial conds `x0` on time grid `tgrid`.
    
    Returns dense output (an array of solutions at each point in `tgrid`) """#
function ode4_dense(F, x0::Number, tgrid)
    solver = GenericRK4(F, x0, tgrid)

    # return a 1D array
    vcat(collect(solver)...)
end
# version for multi-valued functions
function ode4_dense(F, x0::Vector, tgrid)
    solver = GenericRK4(F, x0, tgrid)

    # return a 2D array
    hcat(collect(solver)...)'
end
