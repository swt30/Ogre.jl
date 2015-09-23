# Global configuration settings


"Configuration options"
module config
    # Data directories
    "Location of the data files in this package"
    const datadir = normpath(joinpath(dirname(@__FILE__), "..", "data"))
end

include("config_defaults.jl")
