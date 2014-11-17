module integrator
export BoundaryValues, ODESystem, solve, zero, copy, get_radius
using ogre.common, ogre.constants
import Base.copy, Base.zero

# Mass coordinates and other values

initial_values(vs::ValueSet) = [vs.r, vs.P]
mass_coordinate(vs::ValueSet) = vs.m

# Boundary values and ODE

typealias BoundaryValues ValueSet

type ODESystem
    equations::EquationSet
    boundary_values::BoundaryValues
end

copy(sys::ODESystem) = ODESystem(sys.equations, sys.boundary_values)

type ODESolution
    t::Vector{Float64}
    y::Array{Float64, 2}
end

zero(::Type{ODESolution}) = ODESolution([0.], [0. 0.])

# Stepper function for doing a single step of the integration
function solve_step(ode::ODESystem, t_grid::Vector{Float64})
    ode_func(t::Float64, y::Vector{Float64}) = call(ode.equations, t, y)

    boundary = ode.boundary_values
    y_start::Vector{Float64} = initial_values(boundary)
    t_start::Float64 = mass_coordinate(boundary)

    yn = y_start
    tn = t_start

    produce((t_start, y_start))
    for tnext in t_grid[2:end]
        yn = ode4(ode_func, yn, tn, tnext)
        tn = tnext
        if should_halt(tn, yn)
            # give us one more value and then stop
            produce((tn, yn))
            break
        else
            produce((tn, yn))
        end
    end
end

function should_halt(t, y)
    # if the radius drops below zero, we halt there
    if y[1] < 0
        return true
    else
        return false
    end
end

# Functions to actually solve an ODESystem

function solve_prep(ode::ODESystem, t_grid::Vector{Float64})
    # array setup
    n_points = length(t_grid)
    t = fill(NaN, n_points)
    y = fill(NaN, (n_points, 2))
    solution = zero(ODESolution)
    return t, y, solution
end

function solve!(ode::ODESystem, soln::ODESolution, t_grid::Vector{Float64},
                t::Vector{Float64}, y::Array{Float64, 2})
    # make an integrator stepper and run it
    stepper = @task solve_step(ode, t_grid)
    for (i, x) in enumerate(stepper)
        t[i] = x[1]
        y[i, :] = x[2]
    end

    soln.t = t
    soln.y = y

    return soln
end

function solve(ode::ODESystem, t_grid::Vector{Float64})
    t, y, soln = solve_prep(ode, t_grid)
    soln = solve!(ode, soln, t_grid, t, y)
end

# a couple of tools for use in the radius solution
notnan(arr) = ~isnan(arr)
dropna(arr::Array) = filter(notnan, arr)

function solve_for_radius!(system::ODESystem,
                           solution_grid::Vector{Float64},
                           R_bracket::Vector{Float64})
    R_low, R_high = R_bracket
    done = false
    result = zero(ODESolution)

    while !done
        # choose a radius
        R_guess = (R_low + R_high) / 2

        # set up a new planet system with different radius
        system.boundary_values.r = R_guess

        # solve the system
        result = solve(system, solution_grid)
        radii = result.y[:, 1] |> dropna
        final_radius = radii[end]

        if final_radius < 0
            R_low = R_guess
        elseif final_radius > 100
            R_high = R_guess
        else
            done = true
            return R_guess, result
        end
    end
end

function setup_system(M::Float64, R::Float64,
                      P_surface::Float64, structure_equations::EquationSet)
    # boundary conditions and ODE setup
    bv = BoundaryValues(M, R, P_surface)
    system = ODESystem(structure_equations, bv)
end

function get_radius(system::ODESystem,
                    solution_grid::Vector{Float64},
                    R_bracket::Vector{Float64})
    R, _ = solve_for_radius!(system, solution_grid, R_bracket)
    return R::Float64
end

function get_radius(M::Float64,
                    structure_equations::EquationSet,
                    P_surface::Float64,
                    solution_grid::Vector{Float64},
                    R_bracket::Vector{Float64})
    R_guess = mean(R_bracket)
    system = setup_system(M, R_guess, P_surface,
                          structure_equations)
    R::Float64 = get_radius(system, solution_grid, R_bracket)
end

# numerical methods

# a simple RK4 step
function ode4(F::Function, x0::Vector{Float64}, tstart, tend)
    h = tend - tstart
    n_steps = 2
    xnew = copy(x0)

    k = Array(typeof(x0), 4)
    # Beginning derivative
    k[1] = h*F(tstart,        x0)
    # Midstep derivatives
    k[2] = h*F(tstart + h./2, x0 + k[1]./2)
    k[3] = h*F(tstart + h./2, x0 + k[2]./2)
    # Ending derivative
    k[4] = h*F(tstart + h,    x0 + k[3])

    # Integrate
    xnew = x0 + k[1]./6 + k[2]./3 + k[3]./3 + k[4]./4
end


end
