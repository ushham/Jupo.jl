
function lorenz!(du, u, p, t)
    σ, r, b = p
    du[1] = σ * (u[2] - u[1])
    du[2] = - u[1] * u[3] + r * u[1] - u[2]
    du[3] = u[1] * u[2] - b * u[3]
    return nothing
end

function lorenz_j!(du, u, p, t)
    σ, r, b = p
    du[1, :] = [-σ, σ, 0]
    du[2, :] = [(r - u[3]), -1, -u[1]]
    du[3, :] = [u[2], u[1], -b]
    return nothing
end