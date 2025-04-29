
"""
    Function that extends the function `plot` from `Plots`, to plot a UPO.
    # Arguments
    - `ls::CoupledODEs`: As in the `DynamicalSystems` library.
    - `upo::UPO_sol`: UPO solution struct.
    - `axis::Union{Tuple{Int}}`: Projection axes of the UPO to plot.
    - `Δt::Float64`: Timestep for integration, defaults to `0.01`.

    # Returns
    `plot`

"""
function plot(
    sys::CoupledODEs, 
    upo::UPO_sol; 
    axis::Union{Tuple{T, T}, Tuple{T, T, T}}=(1, 2),
    Δt::Float64=0.01,
    kwargs...
    ) where {T<:Integer}
    traj, t = trajectory(sys, upo.period, upo.ic, Δt=Δt)

    # add first point to last point
    push!(traj, traj[1])

    if length(axis) == 3
        plot(traj[:, axis[1]], traj[:, axis[2]], traj[:, axis[3]], kwargs...)
    else
        plot(traj[:, axis[1]], traj[:, axis[2]], kwargs...)
    end
    return plot!()
end


"""
    Function that extends the function `plot` from `Plots`, to plot a UPO.
    # Arguments
    - `ls::CoupledODEs`: As in the `DynamicalSystems` library.
    - `upo::UPO_sol`: UPO solution struct.
    - `axis::Union{Tuple{Int}}`: Projection axes of the UPO to plot.
    - `Δt::Float64`: Timestep for integration, defaults to `0.01`.

    # Returns
    `plot`

"""
function plot(
    sys::System, 
    upo::UPO_sol; 
    axis::Union{Tuple{T, T}, Tuple{T, T, T}}=(1, 2),
    Δt::Float64=0.01,
    kwargs...
    ) where {T<:Integer}

    ls = generate_ode(sys)

    traj, t = trajectory(ls, upo.period, upo.ic, Δt=Δt)

    # add first point to last point
    push!(traj, traj[1])

    if length(axis) == 3
        plot(traj[:, axis[1]], traj[:, axis[2]], traj[:, axis[3]], kwargs...)
    else
        plot(traj[:, axis[1]], traj[:, axis[2]], kwargs...)
    end
    return plot!()
end


"""
    Inplace function to extend the function `plot!` from `Plots`, to plot a UPO.
    # Arguments
    - `ls::CoupledODEs`: As in the `DynamicalSystems` library.
    - `upo::UPO_sol`: UPO solution struct.
    - `axis::Union{Tuple{Int}}`: Projection axes of the UPO to plot.
    - `Δt::Float64`: Timestep for integration, defaults to `0.01`.

    # Returns
    `plot!`
"""
function plot!(
    sys::CoupledODEs, 
    upo::UPO_sol; 
    axis::Union{Tuple{T, T}, Tuple{T, T, T}}=(1, 2),
    Δt::Float64=0.01,
    kwargs...
    )  where {T<:Integer}

    traj, t = trajectory(sys, upo.period, upo.ic, Δt=Δt)

    # add first point to last point
    push!(traj, traj[1])

    if length(axis) == 3
        plot!(traj[:, axis[1]], traj[:, axis[2]], traj[:, axis[3]], kwargs...)
    else
        plot!(traj[:, axis[1]], traj[:, axis[2]], kwargs...)
    end
    return plot!()
end


"""
    Inplace function to extend the function `plot!` from `Plots`, to plot a UPO.
    # Arguments
    - `sys::System`: As in the `DynamicalSystems` library.
    - `upo::UPO_sol`: UPO solution struct.
    - `axis::Union{Tuple{Int}}`: Projection axes of the UPO to plot.
    - `Δt::Float64`: Timestep for integration, defaults to `0.01`.

    # Returns
    `plot!`
"""
function plot!(
    sys::System, 
    upo::UPO_sol; 
    axis::Union{Tuple{T, T}, Tuple{T, T, T}}=(1, 2),
    Δt::Float64=0.01,
    kwargs...
    )  where {T<:Integer}

    ls = generate_ode(sys)

    traj, t = trajectory(ls, upo.period, upo.ic, Δt=Δt)

    # add first point to last point
    push!(traj, traj[1])
    
    if length(axis) == 3
        plot!(traj[:, axis[1]], traj[:, axis[2]], traj[:, axis[3]], kwargs...)
    else
        plot!(traj[:, axis[1]], traj[:, axis[2]], kwargs...)
    end
    return plot!()
end