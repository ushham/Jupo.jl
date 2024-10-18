#=
stabilising_transforms:
- Julia version: 
- Author: ohamilton
- Date: 2024-05-20
=#
include("base.jl")
include("system.jl")

using DynamicalSystems
using LinearAlgebra
using PrettyTables
using NLsolve
using DifferentialEquations


begin
    Base.@kwdef mutable struct Stability_kwargs
        C̃::Float64 = -1.
        α::Float64 = 0.01
        full_loop::Bool = false
    end
end


function construct_subspaces(tds::TangentDynamicalSystem, x_i::Vector{Float64}, kw::Stability_kwargs)
    delta = 0.001
    if kw.full_loop
        _, m = integrate_tgls(tds, x_i)
    else
        _, m = integrate_tgls(tds, x_i[1:end-1], t=delta)
    end

    U, Σ, _ = svd(m)

    num_d = sum(Σ .>= 1)

    unstable_subspace = Array{Float64}(undef, length(x_i)-1, num_d)
    n = 1
    d = 1
    for s in Σ
        if s >= 1
            unstable_subspace[:, d] = U[:, n]
            d += 1
        end
        n += 1
    end
    return transpose(unstable_subspace)
end

function _dx_ds(tds::TangentDynamicalSystem, x_i::Vector{Float64}, kw::Stability_kwargs)
    F_x_val = F_x(tds, x_i)
    
    V_u = construct_subspaces(tds, x_i, kw)

    d, _ = size(V_u)
    I_mat = Matrix(1I, d, d)
    C_ = kw.C̃ .- I_mat
    C_ = transpose(V_u) * C_ * V_u

    C = Matrix(1I, length(x_i)-1, length(x_i)-1) + C_
    return C * F_x_val
end


function _dT_ds(tds::TangentDynamicalSystem, x_i::Vector{Float64}, kw::Stability_kwargs)
    f = dynamic_rule(tds)
    F_x_val, S_x = F_x(tds, x_i, true)

    f_s_x = Array{Float64}(undef, length(x_i)-1)
    f(f_s_x, S_x, current_parameters(tds), 0)

    return - kw.α * dot(f_s_x, F_x_val)
end


function unstable_system!(dx, x, p, t)
    tds = p[1]
    kw = p[2]
    
    dx[1:end-1] = _dx_ds(tds, x, kw)
    dx[end] = _dT_ds(tds, x, kw)
    return nothing
end


function create_subsystem(tds::TangentDynamicalSystem, ic::Vector{Float64}, kw=Stability_kwargs())
    #\\TODO: Need to work out how to use other solver in this section
    return CoupledODEs(unstable_system!, ic, [tds, kw], diffeq=(abstol = 1e-10, reltol = 1e-10))
end


# lorenz_sys = lorenz
# ls = CoupledODEs(lorenz_sys.f, randn(3), [10., 28., 8/3])
# tds = TangentDynamicalSystem(ls, J=lorenz_sys.jac, Q0=Matrix(1.0I, 3, 3))

# kw = Stability_kwargs()


# ic = [-13.76, -19.57, 27., 1.57]
# test = Array{Float64}(undef, 4)

# test[end] = _dT_ds(tds, ic, kw=kw)
# test


# unstable_system!(test, ic, [tds, Stability_kwargs()], 0)
# create_subsystem(tds, ic)