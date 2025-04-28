using Jupo

lorenz_system = lorenz

# Find "figure of 8" UPO, the orbit that does one loop in each wing.
# Through trial and error this orbit has a period of T ≈ 1.55

# The period is added as the final value of the initial condition vector
ic_point = [14.1, 9.5, 38.64]
ic_period = 1.55
ic_on_attractor = [ic_point; ic_period]
upo = find_upo_nm(lorenz_system, ic_on_attractor; print_report=true)

plot(lorenz_system, upo)