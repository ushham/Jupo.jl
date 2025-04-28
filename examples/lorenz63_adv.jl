using Jupo
using DynamicalSystems
using Plots

lorenz_system = lorenz

# If we do not know the initial conditions there is a function that can help find suitable initial conditions, up to a period T≤5.
ics = guess_ic(lorenz; max_period=5.)

upos = []
for ic_try in ics
    push!(upos, find_upo_nm(lorenz_system, ic_try))
end

# To plot the attractor, we create a long transient
ls = generate_ode(lorenz)
traj, t = trajectory(ls, 1000, rand(3), Δt=0.01)
plot(traj[:, 1], traj[:, 2], lc=:lightgrey)

# Plot all upos
for u in upos
    plot!(lorenz, u)
end

plot!()

# For some systems, or initial conditions, finding a UPO can be more tricky
# In these situations there are a number of parameters to tweak

# Number of iterations done (can stop on other conditions)
# Min norm in the correction step before stopping the algorithm
# Maximum error growth between iterations - this allows for increasing the error to move out of a local minimum in the hope of finding a better minimum
# Damping of the step size algorithm
# Maximum number of steps per step size algorithm
# allowing the tensor correction instead of the more simple linear method 

# These are some suggestions to play with, see the docs for more detail

upo_exact = find_upo_nm(
    lorenz, 
    upos[2];
    print_report=true, 
    iterations=100,
    damping_max_steps=200,
    min_norm=1e-16,
    error_flex=10.,
    damping_α=0.5,
    tensor_auto=true
)