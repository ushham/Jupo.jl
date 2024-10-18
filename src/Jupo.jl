module Jupo

include("methods/base.jl")
include("methods/continuation.jl")
include("methods/floquet.jl")
include("methods/initial_conditions.jl")
include("methods/newton_method.jl")
include("methods/plotting.jl")
include("methods/system.jl")
include("methods/tracking.jl")
include("methods/upo.jl")

include("systems/MAOLAM2x2.jl")
include("systems/toy_systems.jl")

end
