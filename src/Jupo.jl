module Jupo

# Base

using DynamicalSystems: 
    TangentDynamicalSystem,
    CoupledODEs,
    reinit!,
    step!,
    current_state,
    trajectory,
    dimension,
    dynamic_rule,
    current_parameters

using LinearAlgebra:
    Adjoint, 
    I,
    eigvals,   
    svd,
    dot,
    ⋅

using Distances:
    norm

using NLsolve:
    nlsolve

using JSON
using PrettyTables:
    pretty_table

import Plots:
    plot,
    plot!


# Import Systems
include("systems/MAOLAM2x2.jl")
include("systems/toy_systems.jl")

# Import methods
include("methods/system.jl")
include("methods/base.jl")
include("methods/continuation.jl")
include("methods/floquet.jl")
include("methods/initial_conditions.jl")
include("methods/newton_method.jl")
include("methods/plotting.jl")

export System, SystemDef
export lorenz, maolam_2x2
export generate_ode, generate_ds, upo_to_dict, dict_to_upo, save_upo, open_upos, UPO_sol
export continue_upo
export floquet_multipliers, floquet_exponents
export filter_ic_guesses, guess_ic
export find_upo_nm
export plot, plot!

end
