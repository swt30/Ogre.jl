# Default values for the integrator


"Default values for the integrator"
module defaults
import Ogre: bar, K, R_earth

"Surface pressure"
const P_surf = 1bar
"Surface temperature"
const T_surf = 300K
"Mass fractions of the different layers"
const mass_fractions = [1.]
"Total number of points to use in the integration"
const total_points = 500
"Radius search bracket"
const R_bracket = [0., 10.] .* R_earth
end
