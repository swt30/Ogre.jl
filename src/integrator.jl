# INTEGRATOR.JL
# Numerical routines for solutions of planetary structure models

# Solutions
#------------------------------------------------------------------------------

"Solve a `PlanetSystem`, writing the solution to a given `PlanetStructure`."
function solve!(sys::PlanetSystem, soln::PlanetStructure)
    # initial values
    x0 = nonmass(sys.boundary_values)
    t0 = mass(sys.boundary_values)
    tgrid = sys.solution_grid

    # solve and write to solution array
    solver = RK4(sys.structure_equations, x0, tgrid)
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

"Do a RK4 step of function `F`, initial conds `x`, step size `dt`"
function ode4_step(F, x::Real, tstart, tend)
    k = Array(eltype(x), 4)
    dt = tend - tstart

    # Beginning derivative
    k[1] = dt*F(tstart,        x)
    # Midstep derivatives
    k[2] = dt*F(tstart + dt/2, x + k[1]/2)
    k[3] = dt*F(tstart + dt/2, x + k[2]/2)
    # Ending derivative
    k[4] = dt*F(tstart + dt,   x + k[3])

    # Integrate
    +(x, k[1]/6, k[2]/3, k[3]/3, k[4]/6)
end
function ode4_step{R<:Real}(F, x::Vector{R}, tstart, tend)
    k = Array(eltype(x), (length(x), 4))
    dt = tend - tstart
    
    # Beginning derivative
    k[:,1] = dt*F(tstart,        x)
    # Midstep derivatives
    k[:,2] = dt*F(tstart + dt/2, x + k[:,1]/2)
    k[:,3] = dt*F(tstart + dt/2, x + k[:,2]/2)
    # Ending derivative
    k[:,4] = dt*F(tstart + dt,   x + k[:,3])

    # Integrate
    +(x, k[:,1]/6, k[:,2]/3, k[:,3]/3, k[:,4]/6)
end


# types to handle solving
abstract IntegratorMethod
abstract FixedStepIntegrator <: IntegratorMethod
"RK4 solver"
abstract RK4Integrator <: FixedStepIntegrator

abstract IntegratorState
immutable RK4IntegratorState{T} <: IntegratorState
    tindex::Int
    x::T
end

timeindex(s::IntegratorState) = s.tindex
nextindex(s::IntegratorState) = timeindex(s) + 1
x(s::IntegratorState) = s.x
function nextx(s::RK4IntegratorState, solver)
    if timeindex(s) == length(solver.tgrid)
        # can't integrate past the end so leave the x as it is
        x(s)
    else
        # do an integration step
        t = solver.tgrid[timeindex(s)]
        tn = solver.tgrid[nextindex(s)]
        ode4_step(solver.func, x(s), t, tn)
    end
end

"RK4 solver for scalar functions"
immutable ScalarRK4{T<:Real, V<:AbstractVector{Float64}} <: RK4Integrator
    func
    x0::T
    tgrid::V
end
Base.eltype(::Type{ScalarRK4}) = Float64
RK4(func, x0::Real, tgrid) = ScalarRK4(func, Float64(x0), tgrid)

"RK4 solver for vector-valued functions"
immutable VectorRK4{V1<:Vector{Float64}, V2<:AbstractVector{Float64}} <: RK4Integrator
    func
    x0::V1
    tgrid::V2
end
Base.eltype(::Type{VectorRK4}) = Vector{Float64}
RK4(func, x0::AbstractVector, tgrid) = VectorRK4(func, Float64[x0...], tgrid)

# iterator setup
Base.start(solver::IntegratorMethod) = RK4IntegratorState(1, solver.x0)

# step function
function Base.next(solver::RK4Integrator, state)
    (state.x, RK4IntegratorState(nextindex(state), nextx(state, solver)))
end

# termination condition
function Base.done(solver::FixedStepIntegrator, state)
    timeindex(state) > length(solver.tgrid)
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
    
    Returns dense output (an array of solutions at each point in `tgrid`) """
function ode4_dense(F, x0::Real, tgrid)
    solver = RK4(F, x0, tgrid)

    # return a 1D array
    collect(solver)
end
# version for multi-valued functions
function ode4_dense(F, x0::Vector, tgrid)
    solver = RK4(F, x0, tgrid)

    # return a 2D array
    hcat(collect(solver)...)'
end
