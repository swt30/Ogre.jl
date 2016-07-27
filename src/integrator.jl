# Numerical solutions of planetary structure models


# Solver types

abstract Integrator
abstract FixedStepIntegrator <: Integrator

typealias FloatOrVecFloat Union{Float64, Vector{Float64}}

"RK4 solver"
type RK4{ResultType<:FloatOrVecFloat} <: FixedStepIntegrator
    func
    y0::ResultType
    tgrid::AbstractVector
    ks::Vector{ResultType}
    y::ResultType
end
RK4(f, y0::FloatOrVecFloat, tgrid) = RK4(f, y0, tgrid, [zero(y0) for _=1:4], copy(y0))
RK4(f, y0::Real, tgrid) = RK4(f, Float64(y0), tgrid)
RK4{R<:Real}(f, y0::Vector{R}, tgrid) = RK4(f, collect(Float64, y0), tgrid)

typealias ScalarRK4 RK4{Float64}
typealias VectorRK4 RK4{Vector{Float64}}


# Working with integrator and state types

function reset!(integrator::ScalarRK4)
    integrator.y = integrator.y0
    integrator.ks = [zero(integrator.y0) for _=1:4]
    return nothing
end
function reset!(integrator::VectorRK4)
    integrator.y[:] = integrator.y0
    integrator.ks = [zero(integrator.y0) for _=1:4]
    return nothing
end
nvars(integrator) = length(integrator.y0)
ntimes(integrator) = length(integrator.tgrid)
Base.length(integrator::RK4) = ntimes(integrator)
time(integrator, n) = integrator.tgrid[n]
nexttime(integrator, n) = integrator.tgrid[n + 1]
y(integrator::Integrator) = integrator.y
y0(integrator) = integrator.y0
atfinaltime(integrator, n) = (n == ntimes(integrator))
pastfinaltime(integrator, n) = (n > ntimes(integrator))
Base.eltype(::Type{ScalarRK4}) = Float64
Base.eltype(::Type{VectorRK4}) = Vector{Float64}


# Numerical steps

const rk4_dy_coeffs = [1/3, 1/6, 1/6, 1/3]

"Do a RK4 step of function `F`, initial conds `x0` and `t0`, step size `dt`"
function ode4_step!(i::ScalarRK4, t0::Float64, t1::Float64)
    let ks = i.ks, F = i.func, order = 4, dt = t1 - t0

        ks[1] = dt*F(t0,        i.y)
        ks[2] = dt*F(t0 + dt/2, i.y + ks[1]/2)
        ks[3] = dt*F(t0 + dt/2, i.y + ks[2]/2)
        ks[4] = dt*F(t0 + dt,   i.y + ks[3])

        for d = 1:order
            i.y += ks[d] * rk4_dy_coeffs[d]
        end

        nothing
    end
end
function ode4_step!(i::VectorRK4, t0::Float64, t1::Float64)
    let ks = i.ks, y = i.y, F = i.func, n = nvars(i), order = 4, dt = t1 - t0

        ks[1] = dt*F(t0,        y)
        ks[2] = dt*F(t0 + dt/2, y + ks[1]/2)
        ks[3] = dt*F(t0 + dt/2, y + ks[2]/2)
        ks[4] = dt*F(t0 + dt,   y + ks[3])

        for d = 1:order, j=1:n
            y[j] += ks[d][j] * rk4_dy_coeffs[d]
        end

        nothing
    end
end
function ode4_step!(i::Integrator, n)
    let t = time(i, n), tn = nexttime(i, n)
        ode4_step!(i, t, tn)
    end

    nothing
end
function nexty!(i::Integrator, n)
    if atfinaltime(i, n)
        # can't integrate past the end point
        nothing
    else
        # do an integration
        ode4_step!(i, n)
    end
end

# iterator setup
Base.start(rk4::RK4) = (reset!(rk4); return 1)
# stepper
function Base.next(rk4::RK4, n)
    yout = copy(y(rk4))
    nexty!(rk4, n)
    return yout, n + 1
end
# termination condition
Base.done(fs::FixedStepIntegrator, n) = pastfinaltime(fs, n)

# for planets, if a value becomes unphysical (such as r<0 or P<0) then the
# structure equations will yield zeroes for each remaining step. This allows us
# to step "over" the core and leaves us with a slightly negative value at the
# centre. This negative sign signals the radius search function that we should
# increase our radius.


# Dense solvers

""" Solve the ODE dx/dt = `F` with initial conds `x0` on time grid `tgrid`.

    Returns dense output (an array of solutions at each point in `tgrid`) """
function ode4_dense(F, y0::Real, tgrid)
    integrator = RK4(F, y0, tgrid)

    # return a 1D array
    [y for y in integrator]
end
# version for multi-valued functions
function ode4_dense(F, x0::Vector, tgrid)
    integrator = RK4(F, x0, tgrid)

    # return a 2D array
    [y[i] for y in integrator, i in 1:nvars(integrator)]
end


# Solutions

"Solve a `PlanetSystem`, writing the solution to a given `PlanetStructure`."
function solve!(sys::PlanetSystem, soln::PlanetStructure)
    # initial values
    y0 = nonmass(sys.boundary_values)
    tgrid = sys.solution_grid

    # solve and write to solution array
    solver = RK4(sys.structure_equations, y0, tgrid)
    for (i, y) in enumerate(solver)
        setmass!(soln, i, tgrid[i])
        setnonmass!(soln, i, y)
    end

    return soln
end

"Generate a new `PlanetStructure` object and solve a `PlanetSystem`"
function solve(sys::PlanetSystem)
    soln = blank_structure(sys)
    solve!(sys, soln)

    return soln
end

# Radius searching

# helper functions
"Has the radius passed the central point?"
hit_the_centre(R::Radius) = R < 0
hit_the_centre(ps::PlanetStructure) = hit_the_centre(radius(centre(ps)))
"Is the radius not yet close enough to the centre?"
not_far_enough(R::Radius) = R > 100
not_far_enough(ps::PlanetStructure) = not_far_enough(radius(centre(ps)))
unacceptable(ps) = hit_the_centre(ps) || not_far_enough(ps)
"Is the structural solution acceptable?"
acceptable(system, structure) = !unacceptable(structure)

"Update the radius search bracket using a bisection search"
function refine_r_bracket!(system::PlanetSystem, previousresult::PlanetStructure)
    if hit_the_centre(previousresult)
        minR = currentradiusguess(system)
        maxR = maxradius(system)
        refine_r_bracket!(system, minR, maxR)
    elseif not_far_enough(previousresult)
        minR = minradius(system)
        maxR = currentradiusguess(system)
        refine_r_bracket!(system, minR, maxR)
    end
end

"Update the planet's boundary conditions using the results of a test solution"
function refine_boundary_conditions!(system::PlanetSystem, previousresult::PlanetStructure)
    refine_r_bracket!(system, previousresult)
    refine_boundary_r!(system)
    refine_surface!(system)
end

import Plots: plot
Plots.pyplot(show=true)
""" Repeatedly solve a `PlanetSystem` until its radius is suitable

    Recursively search the radius search bracket of a `PlanetSystem` to
    get a result that has an acceptable error at the centre. Returns a
    `PlanetSystem` with the appropriate radius. """
function converge!(system)
    # first make sure the system has the correct boundary conditions to start with
    refine_boundary_r!(system)
    refine_surface!(system)

    # then get the structure for this radius guess
    result = solve(system)
    plot(result)

    # and recursively iterate
    converge!(system, result)
    return system
end
function converge!(system, result)
    # we have already done an iteration and have a result to work with
    if !acceptable(system, result)
        # if diff(system.radius_search_bracket)[1] < 1 && radius(centre(result)) > 0
        # info("Did not converge fully")
        # return system
        # end

        # need to refine and try again
        refine_boundary_conditions!(system, result)
        solve!(system, result)
        plot(result)
        converge!(system, result)
    else
        # we're done
        return system
    end
end

# basic functions for geting radii and structures
"Get the correct internal structure for a `PlanetSystem`"
find_structure!(system) = solve(converge!(system))
"Find the radius of a `PlanetSystem` by solving the system"
find_radius!(system) = radius(surface(find_structure!(system)))

"Solve a `PlanetSystem` and give both the total radius and the structure"
function find_structure_and_radius!(system::PlanetSystem)
    struct = find_structure!(system)
    r = radius(surface(struct))

    struct, r
end

# higher level functions for doing MR diagrams
"Create and solve for a `PlanetSystem`'s radius'"
function R(M::Mass, structure_equations::EquationSet, Psurf::Pressure,
           solution_grid, R_bracket)

    R_guess = mean(R_bracket)
    system = PlanetSystem(M, R_guess, Psurf, structure,
                          solution_grid, R_bracket)

    radius = find_radius!(system)
end

"Find the radius of a solid sphere of mass `M` using an given `EOS`."
function R(M, eos; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1.
    Rscale = in_earth_units ? 1/R_earth : 1.

    sys = DefaultPlanetSystem(M * Mscale, eos)
    radius = find_radius!(sys) * Rscale
end

# temperature-dependent version
function R(M, eos, Cₚ; in_earth_units=false)
    Mscale = in_earth_units ? M_earth : 1.
    Rscale = in_earth_units ? 1/R_earth : 1.

    sys = DefaultPlanetSystem(M * Mscale, eos, Cₚ)
    radius = find_radius!(sys) * Rscale
end

# vectorized form of the above
function R{T<:Real}(ms::AbstractVector{T}, args...; in_earth_units=false)
    R_withEOS(M) = R(M, args...; in_earth_units=in_earth_units)
    radii = map(R_withEOS, ms)
end
