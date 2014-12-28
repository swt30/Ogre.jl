module integrator
export BoundaryValues, PlanetSystem, PlanetStructure
export zero, copy, find_radius, find_radius!, solve
using ogre.common, ogre.constants
import Base.copy, Base.zero

#= Mass coordinates and other values =#

initial_values{T<:Real}(vs::ValueSet{T}) = [vs.r, vs.P]::Vector{T}
mass_coordinate(vs::ValueSet) = vs.m::Real

#= Types that define the planetary structure and the problem =#

typealias BoundaryValues ValueSet

type PlanetSystem{T<:Real}
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

zero(::Type{PlanetStructure}) = PlanetStructure([0.], [0. 0.])

#= Functions to actually solve the structure, given boundary conditions =#

# stepper task for doing the integration
function stepper(sys::PlanetSystem)
    # should be able to overload the call method here eventually
    ode_func{T<:Real}(t::T, y::Vector{T}) = callfunc(sys.structure_equations,
                                                     t, y)

    boundary = sys.boundary_values
    ystart = initial_values(boundary)
    tstart = mass_coordinate(boundary)
    tn = sys.solution_grid[1]

    F, yn, k = ode4_setup(ode_func, ystart)

    produce((tstart, yn))
    for tnext in sys.solution_grid[2:end]
        ode4!(F, yn, k, tn, tnext)
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

# callback function to tell if the integration should stop
function should_halt{T<:Real}(t::T, y::Vector{T})
    # if the radius or mass drops below zero, we halt there
    if y[1] < 0 || t < 0
        return true
    else
        return false
    end
end

# generating blank solutions
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
    # make an integrator stepper and run it
    steppertask = @task stepper(sys)
    for (i, x) in enumerate(steppertask)
        soln.m[i] = x[1]
        soln.y[i, :] = x[2]
    end
end

# generate a new solution object and solve
function solve(sys::PlanetSystem)
    soln = blank_structure(sys)
    solve!(sys, soln)

    soln
end

#= functions to get a radius for a structure =#

# helper functions
notnan(arr) = ~isnan(arr)
dropnans(arr::Array) = filter(notnan, arr)

# iterate to get a radius that is suitable
function find_radius!(system::PlanetSystem)
    R_low, R_high = system.radius_search_bracket
    done = false
    result = blank_structure(system)

    while !done
        # choose a radius
        R_guess = (R_low + R_high) / 2
        system.boundary_values.r = R_guess

        # solve the system for that radius
        result = solve(system)
        radii = result.y[:, 1] |> dropnans
        final_radius = radii[end]

        # rinse and repeat
        if final_radius < 0
            R_low = R_guess
        elseif final_radius > 100
            R_high = R_guess
        else
            done = true
            return R_guess
        end
    end
end

# set up a planetary system
function setup_find_radius{T<:Real}(M::T,
                                    R::T,
                                    P_surface::T,
                                    struct::EquationSet,
                                    solution_grid::Vector{T},
                                    R_bracket::Vector{T})
    # boundary conditions and ODE setup
    bv = BoundaryValues(M, R, P_surface)
    system = PlanetSystem(M, struct, bv, solution_grid, R_bracket)
end

# create a system and find an appropriate radius
function find_radius{T<:Real}(M::T,
                              structure::EquationSet,
                              P_surface::T,
                              solution_grid::Vector{T},
                              R_bracket::Vector{T})
    R_guess = mean(R_bracket)
    system = setup_system(M, R_guess, P_surface, structure,
                          solution_grid, R_bracket)
    radius = zero(PlanetStructure)
    R::T = find_radius!(system)
end

#= NUMERICAL METHODS =#

# set up arrays for fast ODE solving
function ode4_setup{T<:Real}(F::Function, x::Vector{T})
    xnew::Vector{T} = copy(x)
    k::Matrix{T} = zeros(length(xnew), 4)
    F, xnew, k
end

# a single RK4 step
function ode4!{T<:Real}(F::Function, x::Vector{T}, k::Matrix{T},
                        tstart::T, tend::T)
    h::T = tend - tstart

    # Beginning derivative
    k[:,1] = h*F(tstart,        x)
    # Midstep derivatives
    k[:,2] = h*F(tstart + h./2, x + k[:,1]./2)
    k[:,3] = h*F(tstart + h./2, x + k[:,2]./2)
    # Ending derivative
    k[:,4] = h*F(tstart + h,    x + k[:,3])

    # Integrate and change the x array
    x[:] += k[:,1]./6 + k[:,2]./3 + k[:,3]./3 + k[:,4]./4
end

end
