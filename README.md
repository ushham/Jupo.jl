# Jupo

[![Build Status](https://github.com/ushham/Jupo.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ushham/Jupo.jl/actions/workflows/CI.yml?query=branch%3Amain)

A julia package for numerically finding Unstable Periodic Orbits (UPOs) in dynamical systems.

This library accompanies the study currently under review:

Hamilton, O., Demaeyer, J., Crucifix, M. & Vannitsem, S. (2025) Using Unstable Periodic Orbits to Understand Blocking Behaviour in a Low Order Land-Atmosphere Model [under review] [https://arxiv.org/abs/2503.02808](https://arxiv.org/abs/2503.02808)

It is planned to add this library to [ChaosTools.jl](https://github.com/JuliaDynamics/ChaosTools.jl), see issue: 

### Usage Examples
Given a dynamical system, such as the Lorenz 63 system:
$$
\begin{aligned}
& \frac{\mathrm{d} x}{\mathrm{~d} t}=\sigma(y-x), \\
& \frac{\mathrm{d} y}{\mathrm{~d} t}=x(\rho-z)-y, \\
& \frac{\mathrm{~d} z}{\mathrm{~d} t}=x y-\beta z .
\end{aligned}
$$ 

and its Jacobian:
$$
J=\begin{bmatrix}
-\sigma & \sigma & 0 \\
\rho - z & -1 & -x \\
y & x & -\beta
\end{bmatrix}
$$

We define a `System`, as the collection of the dynamical rule, the jacobian, and some other information about the system:
```julia
lorenz = System(
        name="Lorenz", 
        f=lorenz!, 
        jac=lorenz_j!, 
        ndim=3, 
        p=[10., 28., 8/3],
        p_names=["sig", "r", "b"]
        )
end
```

this is passed to the method of choice for finding the UPOs, along with a `Vector` containing the initial conditions and the initial period: `ic = [x_0; T]`. The current methods availible are:
- Newton method

Methods to be added:
- Stabilising Transforms

```julia
upo = find_upo_nm(lorenz, ic)
```

This will produce a `UPO_sol` struct, which contains:
- success (whether the algorithm ran successfully)
- ic (initial conditions)
- period
- system (summary of system parameters)
- floquet_multipliers

An example for the Lorenz system, for XXX UPOs is shown below: