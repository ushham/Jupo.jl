"""
    Functions for tracking which UPOs a trajectory is close to.
"""

include("system.jl")
include("base.jl")

using DynamicalSystems
using LinearAlgebra
using PrettyTables
using NLsolve


function format_upos(ls::CoupledODEs, upos::Vector{Any})::Vector{Matrix{Float64}}
    """
        Function generates the trajectory for each UPO input
        Returns a list of trajectories
    """
    vec_out = Vector{Matrix{Float64}}(undef, length(upos))
    for (i, u) in enumerate(upos)
        traj, _ = trajectory(ls, u.period, u.ic, Δt=0.1)
        vec_out[i] = Matrix(traj)
    end
    return vec_out
end


function norms_from_point(upo_traj::Vector{Matrix{Float64}}, point::SVector)::Vector{Float64}
    norms = Vector{Float64}(undef, length(upo_traj))
    for (i, m) in enumerate(upo_traj)
        norms[i] = minimum(norm.(eachrow(m .- point')))
    end
    return norms
end


function closest_upo(ls::CoupledODEs, upos::Vector{Any}, traj::StateSpaceSet)

    upo_mat = format_upos(ls, upos)

    upos_ix = Vector{Int64}(undef, length(traj))
    for (i, p) in enumerate(traj)
        norms = norms_from_point(upo_mat, p)
        upos_ix[i] = argmin(norms)
    end
    return upos_ix
end


