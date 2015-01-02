module Ogre

# add the package directory to the path so that we can find submodules
SRC_PATH = Base.source_path() |> dirname
if !(SRC_PATH in LOAD_PATH)
    push!(LOAD_PATH, Base.source_path() |> dirname)
end

# bring in submodules
include("common.jl")
include("constants.jl")
include("eos.jl")
include("structure.jl")
include("integrator.jl")

using Reexport
@reexport using .Common, .Constants, .Eos, .Structure, .Integrator

end
