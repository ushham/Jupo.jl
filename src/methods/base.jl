include("system.jl")
include("floquet.jl")

using DynamicalSystems
using Distances
import JSON


"""
    A `UPO_sol` struct that represents a periodic orbit (PO).

    # Fields
    - `success::Bool`: Whether the generating algorithm succeeded.
    - `ic::Vector{AbstractFloat}`: The initial condition of the PO.
    - `period::AbstractFloat`: The period of the PO.
    - `error::AbstractFloat`: The error between the initial condition and the final point of the PO
    - `system::SystemDef`: System definition of the PO
    - `floquet_multipliers::Union{Vector{ComplexF64}, Nothing}`: Floquet multipliers of the PO.
"""
Base.@kwdef struct UPO_sol{T<:AbstractFloat}
    success::Bool
    ic::Vector{T}
    period::T
    error::T
    system::SystemDef
    floquet_multipliers::Union{Vector{Complex{T}}, Vector{T}, Nothing} = nothing
end


"""
    A `Floquet_exp` struct that represents the real and imaginary components of the Floquet exponents of a periodic orbit: `Λ_j=exp(T (μ_j + i ω_j))`.

    # Fields
    - `μ::Vector{AbstractFloat}`: The real component of the exponent.
    - `ω::Vector{AbstractFloat}`: The imaginary component of the exponent.

"""
struct Floquet_exp{T<:AbstractFloat}
    μ::Vector{T}
    ω::Vector{T}
end


"""
    Takes a UPO solution and creates a vector containing initial condition in space, with the UPO period: `v = (x_0, T)`.

    # Arguments
    - `sol::UPO_sol`.

    # Returns
    - `Vector{AbstractFloat}`: A vector of the initial condition with the period appended to the end.

"""
function upo_sol_to_sead(sol::UPO_sol)
    output = [sol.ic; sol.period]
    return output
end


"""
    Integrates the tangent linear system in parallel with an initial condition.

    # Arguments
    - `tds::TangentDynamicalSystem`: As in the DynamicalSystems library
    - `x_i::Vector{AbstractFloat}`: Initial condition to begin integration.
    - `t::Union{AbstractFloat, Nothing}=nothing`: Integration time. If `nothing`, this assumes the last element of the vector `x_i` is the time to integrate over.
    - `basis::Union{Matrix{AbstractFloat}, Nothing}=nothing`: Initial condition for the perturbations of the tangent liner system. If nothing it reverts to the identity matrix.
    - `return_x::Bool=true`: Switch to determin if the perturbations from the tangent linear system are returned or not, defaults to `true`.

    # Returns
    If `return_x` is `true`:
    - `Tuple{Vector{AbstractFloat}, Matrix{AbstractFloat}}`: The resulting point in space after the integration and the matrix containing the resulting perturbations of the tangent linear system resulting from the integration.
    
    If `return_x` is false:
    - `Vector{AbstractFloat}`: The resulting point in space after the integration.
"""
function integrate_tgls(
        tds::TangentDynamicalSystem,
        x_i::Vector{T};
        t::Union{T, Nothing}=nothing,
        basis::Union{Matrix{T}, Nothing}=nothing,
        return_x::Bool=true
    ) where {T<:AbstractFloat}

    ndim = dimension(tds)
    if isnothing(t)
        t = x_i[end]
        x_i = x_i[1:ndim]
    end

    if t < 0
        s_x, m = x_i, Matrix(1.0I, length(x_i), length(x_i))
    else
        if isnothing(basis)
            reinit!(tds, x_i)
        else
            reinit!(tds, x_i, Q0=basis)
        end

        step!(tds, t, true)
    end

    if return_x
        return current_state(tds.ds)[:, 1], current_state(tds.ds)[:, 2:end] 
    else
        return current_state(tds.ds)[:, 2:end]
    end
end


"""
    In place integration of the tangent linear system in parallel with an initial condition.

    # Arguments
    - `tds::TangentDynamicalSystem`: As in the DynamicalSystems library
    - `x_i::Vector{AbstractFloat}`: Initial condition to begin integration.
    - `modified_x::Vector{AbstractFloat}`: Vector to be modified by the function.
    - `t::Union{AbstractFloat, Nothing}`: Integration time.
    - `basis::Union{Matrix{AbstractFloat}, Nothing}`: Initial condition for the perturbations of the tangent liner system.

    # Returns
    `nothing`
"""
function integrate_tgls!(
        tds::TangentDynamicalSystem,
        x_i::Vector{T},
        modified_x::Vector{T},
        t::T,
        basis::Matrix{T}
    ) where {T<:AbstractFloat}

    reinit!(tds, x_i, Q0=basis)
    step!(tds, t, true)
    nx = size(current_state(tds.ds), 2)
    modified_x .= view(current_state(tds.ds), :, 1)
    basis .= view(current_state(tds.ds), :, 2:nx)
end


"""
    In place integration of the dynamical system.

    # Arguments
    - `ds::CoupledODEs`: As in the DynamicalSystems library
    - `x_i::Vector{AbstractFloat}`: Initial condition to begin integration.
    - `modified_x::Vector{AbstractFloat}`: Vector to be modified by the function.
    - `t::Union{AbstractFloat, Nothing}`: Integration time.

    # Returns
    `nothing`
"""
function integrate_traj!(
        ds::CoupledODEs,
        x_i::Vector{T},
        modified_x::Vector{T},
        t::T
    ) where {T<:AbstractFloat}
    reinit!(ds, x_i)
    step!(ds, t, true)
    modified_x .= current_state(ds)
end


"""
    Calculates the vector that represents the difference between the initial condition and the resulting point in phase space: `x_0 - F(T, x_0)`

    # Arguments
    - `tds::TangentDynamicalSystem`: As in the DynamicalSystems library
    - `x_i::Vector{AbstractFloat}`: Initial condition to begin integration, with the integration time, or period, stored in the last element of the vector.
    - `return_Sx::Bool=false`: Switch to determine if the function returns the position of the final point in phase space, defulats to `false`.

    # Returns
    If `return_Sx` is `true`:
    - `Tuple{Vector{AbstractFloat}, Vector{AbstractFloat}}`: The difference between the initial condition and the resutling point in phase space after integration, and the final point in phase space.
    
    If `return_x` is false:
    - `Vector{AbstractFloat}`: The difference between the initial condition and the resutling point in phase space after integration.
"""
function F_x(
        tds::Union{TangentDynamicalSystem, CoupledODEs},
        x_i::Vector{T};
        return_Sx::Bool=false
    ) where {T<:AbstractFloat}

    reinit!(tds, x_i[1:end-1])
    step!(tds, x_i[end], true)

    if return_Sx
        return current_state(tds)[:, 1] - x_i[1:end-1], current_state(tds)[:, 1]
    else
        return current_state(tds)[:, 1] - x_i[1:end-1]
    end
end


"""
    Generates a `CoupledODEs` from the `DynamicalSystems` library.

    # Arguments
    - `sys::System`: `struct` containing information on the system properties.
    - `params::Union{Nothing, Vector{AbstractFloat}}=nothing`: Parameters of the system if different from those passed in the `sys` variable. Defulats to `nothing`.
    
    # Returns
    - `CoupledODEs`: System of ODEs as defined in `DynamicalSystems` library.
"""
function generate_ode(
    sys::System; 
    params::Union{Nothing, Vector{T}}=nothing
    )::CoupledODEs where {T<:AbstractFloat}
    if isnothing(params)
        params = sys.p
    end
    ls = CoupledODEs(sys.f, rand(sys.ndim), params, diffeq=(abstol = 1e-16, reltol = 1e-16))
    return ls
end


"""
    Generates a `CoupledODEs` and `TangentDynamicalSystem` from the `DynamicalSystems` library.

    # Arguments
    - `sys::System`: `struct` containing information on the system properties.
    - `params::Union{Nothing, Vector{AbstractFloat}}=nothing`: Parameters of the system if different from those passed in the `sys` variable. Defulats to `nothing`.
    
    # Returns
    - `Tuple{CoupledODEs, TangentDynamicalSystem}`: System of ODEs and the tangent linear system as defined in `DynamicalSystems` library.
"""
function generate_ds(
    sys::System; 
    params::Union{Nothing, Vector{T}}=nothing
    ) where {T<:AbstractFloat}
    ls = generate_ode(sys; params=params)
    tds = TangentDynamicalSystem(ls, J=sys.jac, Q0=Matrix(1.0I, sys.ndim, sys.ndim))

    return ls, tds
end




"""
    Length of a trajectory in phase space.

    # Arguments
    - `traj::Vector{AbstractFloat}`: Trajectory in phase space.
    
    # Returns
    - `AbstractFloat`: Length of the trajectory.
"""
function traj_length(traj)
    arc_len = 0    
    for (t1, t2) in zip(traj[1:end-1], traj[2:end])
        arc_len += norm(t1 - t2)
    end
    return arc_len
end


### Reading/writing outputs ###

"""
    Converts a solution of a periodic solution `UPO_sol` to a dictionary for saving.

    # Arguments
    - `upo::UPO_sol`: Periodic solution.
    
    # Returns
    - `Dict`: Dictionary containing the information of the periodic orbit.
"""
function upo_to_dict(upo::UPO_sol)
    return Dict(
        "success" => upo.success,
        "ic"  => upo.ic,
        "period"  => upo.period, 
        "error"  => upo.error, 
        "floquet_multipliers" => upo.floquet_multipliers,
        "system" => system_to_dict(upo.system),

    )
end


"""
    Converts a dictionary to a `UPO_sol` struct.

    # Arguments
    - `d::Dict`: Dictionary of the key information of the periodic orbit.
    - `system::Union{System, SystemDef, Nothing}=nothing`: Definition of the system of equations of the dynamical systems.
    
    # Returns
    - `UPO_sol`: Periodic orbit solution.
"""
function dict_to_upo(
    d::Dict{String, Any};
    system::Union{System, SystemDef, Nothing}=nothing
    )::UPO_sol

    if isnothing(system)
        system = dict_to_system(d["system"])
    end
    if system isa System
        system = SystemDef(system)
    end
    if haskey(d, "floquet_multipliers")
        fm_dic = d["floquet_multipliers"]
        if eltype(fm_dic) <: AbstractFloat
            fm = map(x -> Float64(x), d["floquet_multipliers"])
        else  
            fm = map(x -> ComplexF64(x["re"] + x["im"] * 1im), d["floquet_multipliers"])
        end 
        return UPO_sol(
            success=d["success"], 
            ic=map(x -> Float64(x), d["ic"]), 
            period=d["period"], 
            error=d["error"], 
            system=system,
            floquet_multipliers=fm
            ) 
    else
        return UPO_sol(
            success=d["success"], 
            ic=map(x -> Float64(x), d["ic"]), 
            period=d["period"], 
            error=d["error"], 
            system=system
            )
    end
end


"""
    Saves a periodic orbit solution `UPO_sol` to a json database.

    # Arguments
    - `upo::Union{UPO_sol, Vector{UPO_sol}}`: Periodic orbit solution, or a vector of periodic orbit solutions.
    - `file_name::Union{String, Nothing, Vector{String}}=nothing`: File name to save the solution to, automatically generates the filename if not given. 
    - `folder::Union{String, Nothing}=nothing`: Folder to save the solution to.
    
    # Returns
    `nothing`
"""
function save_upo(
    upo::Union{UPO_sol, Vector{UPO_sol}}; 
    file_name::Union{String, Nothing, Vector{String}}=nothing, 
    folder::Union{String, Nothing}=nothing
    )
    if upo isa Vector{UPO_sol}
        for (i, u) in enumerate(upo)
            if typeof(file_name) == Vector{String}
                fn = file_name[i]
            else
                fn = nothing
            end
            _write_upo(u, file_name=fn, folder=folder)
        end
    elseif upo isa UPO_sol
        if typeof(file_name) == String
            fn = file_name
        else
            fn = nothing
        end
        _write_upo(upo, file_name=fn, folder=folder)
    end
end


"""
    Hidden function that is used by `save_upo` to write the json file.

    # Arguments
    - `upo::Union{UPO_sol, Vector{UPO_sol}}`: Periodic orbit solution, or a vector of periodic orbit solutions.
    - `file_name::Union{String, Nothing, Vector{String}}=nothing`: File name to save the solution to, automatically generates the filename if not given. 
    - `folder::Union{String, Nothing}=nothing`: Folder to save the solution to.
    
    # Returns
    `nothing`
"""
function _write_upo(
    upo::UPO_sol; 
    file_name::Union{String, Nothing}=nothing, 
    folder::Union{String, Nothing}=nothing
    )
    if isnothing(file_name)
        fn = upo.system.name * "_" * string(round(upo.period, sigdigits=5))
    else
        fn = file_name
    end
    if !isnothing(folder)
        fn = folder * "/" * fn
    end
    if ispath(fn * ".upo")
        i = 0
        not_written = true
        while not_written & (i < 1e2)
            if !ispath(fn * "_" * string(i) * ".upo")
                open(fn * "_" * string(i) * ".upo", "w") do io
                    JSON.print(io, upo_to_dict(upo)) 
                end
                not_written = false
            end
            i += 1
        end
    else
        open(fn * ".upo", "w") do io
            JSON.print(io, upo_to_dict(upo)) 
        end 
    end
end


"""
    Hidden function that is used by `open_upo` to read a saved periodic orbit from a json database.

    # Arguments
    - `loc::String`: File name of the json database.
    - `system::Union{System, SystemDef, Nothing}=nothing`: Information about the dynamical system that the periodic orbit is from, defulats to `nothing`. 

    # Returns
    `UPO_sol`: `struct` containing the periodic orbit information.
"""
function _open_upo(
    loc::String;
    system::Union{System, SystemDef, Nothing}=nothing
    )::Union{UPO_sol, Nothing}

    if split(loc, ".")[end] == "upo"
        loaded_struct = open(loc, "r") do io
            return dict_to_upo(JSON.parse(io), system=system)
        end
    else
        println("Expecting file type upo")
        return nothing
    end
end


"""
    Function to open a saved periodic orbit from a given file.

    # Arguments
    - `loc::String`: File name or folder name containing the files. If a folder name is given all periodic orbits in the file are returned.
    - `min_period::AbstractFloat=0.`: Filters minimum period that is returned if `loc` is the name of a folder of mulitple files. 
    - `max_period::AbstractFloat=1000.`: Filters minimum period that is returned if `loc` is the name of a folder of mulitple files.
    - `system::Union{System, SystemDef, Nothing}=nothing`: Information about the dynamical system that the periodic orbit is from, defulats to `nothing`. 

    # Returns
    `Union{UPO_sol, Vector{UPO_sol}, Nothing}`: Either an array of, or a single `struct` containing the periodic orbit information.
"""
function open_upos(
    loc::String; 
    min_period::T=0., 
    max_period::T=1000.,
    system::Union{System, SystemDef, Nothing}=nothing
    )::Union{UPO_sol, Vector{UPO_sol}, Nothing} where {T<:AbstractFloat}

    upos = nothing
    if isfile(loc)
        temp_upos = _open_upo(loc, system=system)
        if !isnothing(temp_upos)
            if (temp_upos.period >= min_period) & (temp_upos.period <= max_period)
                upos = temp_upos
            end
        end
    end
    if isdir(loc)
        upos = Vector{UPO_sol}(undef, 0)
        files = readdir(loc)
        for f in files
            temp_upos = _open_upo(loc * "/" * f, system=system)
            if !isnothing(temp_upos)
                if (temp_upos.period >= min_period) & (temp_upos.period <= max_period)
                    push!(upos, temp_upos)
                end
            end
        end
    end
    return upos
end
