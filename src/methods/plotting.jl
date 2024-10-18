
function plot(
    ls::CoupledODEs, 
    upo::UPO_sol; 
    axis::Union{Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64}}=(1, 2),
    Δt::Float64=0.01
    )
    traj, t = trajectory(ls, upo.period, upo.ic, Δt=Δt)
    if length(axis) == 3
        plot(traj[:, axis[1]], traj[:, axis[2]], traj[:, axis[3]])
    else
        plot(traj[:, axis[1]], traj[:, axis[2]])
    end
    return plot!()
end

function plot!(
    ls::CoupledODEs, 
    upo::UPO_sol; 
    axis::Union{Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64}}=(1, 2),
    Δt::Float64=0.01
    )

    traj, t = trajectory(ls, upo.period, upo.ic, Δt=Δt)
    if length(axis) == 3
        plot!(traj[:, axis[1]], traj[:, axis[2]], traj[:, axis[3]])
    else
        plot!(traj[:, axis[1]], traj[:, axis[2]])
    end
    return plot!()
end