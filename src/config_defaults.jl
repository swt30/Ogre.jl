# Default values for the integrator


"Default values for the integrator"
module defaults
    import Ogre: R_earth

    "Surface pressure"
    const P_surf = 1e5
    "Surface temperature"
    const T_surf = 300
    "Mass fractions of the different layers"
    const mass_fractions = [1.]
    "Total number of points to use in the integration"
    const total_points = 200
    "Radius search bracket"
    const R_bracket = [0., 10.] .* R_earth
end
