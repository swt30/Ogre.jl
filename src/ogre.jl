module ogre

# add the package directory to the path so that we can find the
SRC_PATH = Base.source_path() |> dirname
if !(SRC_PATH in LOAD_PATH)
    push!(LOAD_PATH, Base.source_path() |> dirname)
end

# bring in submodules
include("common.jl")
include("constants.jl")
include("eos.jl")
include("integrator.jl")
include("structure.jl")

# bring names into the namespace and then reexport them
using Reexport
@reexport using ogre.common
@reexport using ogre.constants
@reexport using ogre.eos
@reexport using ogre.integrator
@reexport using ogre.structure

end
