include("system.jl")
include("base.jl")

using DynamicalSystems
using LinearAlgebra
using PrettyTables
using NLsolve


begin

    function objective_function(x::Vector{T})::T where {T<:AbstractFloat}
        return norm(x, 2)^2
    end

    function grad_obj_function(x::Vector{T})::Vector{T} where {T<:AbstractFloat}
        return 2*x
    end

end


function construct_bounds(arr, tup)
    for (ix, val) in tup
        arr[ix] = val
    end
    return arr
end



"""
    Hidden function that allocates variables that are re-used in each loop of the algorithm.

    # Fields
    - `ndim:Int`: Dimension of the dyamical system.
    
    # Returns
    `Tuple{
        
    }`
"""
function _prealocate_variables(ndim)
    s_x = zeros(ndim)
    lhs_mat = zeros(ndim+1, ndim+1)
    rhs_vec = zeros(ndim+1)
    f_x = Array{Float64}(undef, ndim)
    f_s_x = Array{Float64}(undef, ndim)
    return s_x, lhs_mat, rhs_vec, f_x, f_s_x
end


"""
    Hidden function that generates a `TangentDynamicalSystem` that contains only one basis vector.

    # Arguments
    - `ls::CoupledODEs`: As in the `DynamicalSystems` library.
    - `tds::TangentDynamicalSystem`: As in the `DynamicalSystems` library.
    - `Q0::Union{Nothing, Matrix}=nothing`: Determines the basis of perturbations for the tangent linear system. If `nothing`, a random vector of length one is generated.

    # Returns
    `TangentDynamicalSystem`
"""
function _construct_tds_one_vec(
    ls::CoupledODEs, 
    tds::TangentDynamicalSystem; 
    Q0::Union{Nothing, Matrix}=nothing
    )::TangentDynamicalSystem

    if isnothing(Q0)
        Q0 = rand(dimension(ls), 1)
    end
    return TangentDynamicalSystem(ls, J=tds.J, Q0=Q0)
end


"""
    Inplace function that updates the matrix `A` and the vector `g` in `A Δx = g`, the linearised UPO problem. 

    # Arguments
    - `tds::TangentDynamicalSystem`: As in the `DynamicalSystems` library, tangent linear system of the model.
    - `x_i::Vector{AbstractFloat}`: Initial condition or initial guess for finding the UPO.
    - `t::AbstractFloat`: Initial guess of the period of the uPO.
    - `ndim::Int`: Dimension of the system.
    - `s_x::Vector{AbstractFloat}`: Vector of the difference between the initial guess `x_i`, and the point in the phase space after integration time of `t` from the initial point `x_i`.
    - `lhs_mat::Matrix{AbstractFloat}`: Matrix `A` that describes the linearised UPO problem. This variable is over written so any matrix of the correct dimension can be passed.
    - `rhs_vec::Vector{AbstractFloat}`: Vector `g` that describes the linearised UPO problem. This variable is over written so any vector of the correct dimension can be passed.
    - `f_x::Array{AbstractFloat}`: The tendencies evaluated at the initial point `x_i`, this variable is over written so any vector of the correct dimension can be passed.
    - `f_s_x::Array{AbstractFloat}`: The tendencies evaluated at the point in phase space defined by the initial point `x_i` integrated for time `t`, this variable is over written so any vector of the correct dimension can be passed.
    - `f::Function`: Tendencies of the model `dx/dt=f(x)`.

    # Returns
    `nothing`

"""
function construct_matrix!(;
        tds::TangentDynamicalSystem, 
        x_i::Vector{T},
        t::T,
        ndim::Int,
        s_x::Vector{T},
        lhs_mat::Matrix{T},
        rhs_vec::Vector{T},
        f_x::Array{T},
        f_s_x::Array{T},
        f::Function
    ) where {T<:AbstractFloat}
    m = Matrix(1.0I, ndim, ndim)
    integrate_tgls!(tds, x_i, s_x, t, m)
    
    f(f_s_x, s_x, current_parameters(tds), t)
    f(f_x, x_i, current_parameters(tds), 0)
    

    lhs_mat[1:ndim, 1:ndim] .= m - I
    lhs_mat[1:ndim, end] .= f_s_x
    lhs_mat[end, 1:ndim] .= f_x

    @. rhs_vec[1:ndim] = x_i - s_x
    rhs_vec[end] = 0
    return nothing
end


"""
    Used singular value decomposition to find the solution to`A Δx = g`, the linearised UPO problem. 

    # Arguments
    - `tds::TangentDynamicalSystem`: As in the `DynamicalSystems` library, tangent linear system of the model.
    - `x_i::Vector{AbstractFloat}`: Initial condition or initial guess for finding the UPO.
    - `t::AbstractFloat`: Initial guess of the period of the uPO.
    - `ndim::Int`: Dimension of the system.
    - `s_x::Vector{AbstractFloat}`: Vector of the difference between the initial guess `x_i`, and the point in the phase space after integration time of `t` from the initial point `x_i`.
    - `lhs_mat::Matrix{AbstractFloat}`: Matrix `A` that describes the linearised UPO problem. This variable is over written so any matrix of the correct dimension can be passed.
    - `rhs_vec::Vector{AbstractFloat}`: Vector `g` that describes the linearised UPO problem. This variable is over written so any vector of the correct dimension can be passed.
    - `f_x::Array{AbstractFloat}`: The tendencies evaluated at the initial point `x_i`, this variable is over written so any vector of the correct dimension can be passed.
    - `f_s_x::Array{AbstractFloat}`: The tendencies evaluated at the point in phase space defined by the initial point `x_i` integrated for time `t`, this variable is over written so any vector of the correct dimension can be passed.
    - `f::Function`: Tendencies of the model `dx/dt=f(x)`.

    # Returns
    `Tuple{
        Vector{Float64}, 
        Vector{Float64}
    }`: 
    Given `U Σ V* Δx= A Δx=g`, the inverted linear problem `V(U*g) / Σ`, with the singular values Σ

"""
function svd_problem(;
        tds::TangentDynamicalSystem, 
        x_i::Vector{T},
        t::T,
        ndim::Int,
        s_x::Vector{T},
        lhs_mat::Matrix{T},
        rhs_vec::Vector{T},
        f_x::Array{T},
        f_s_x::Array{T},
        f::Function
    ) where {T<:AbstractFloat}

    construct_matrix!(
        tds=tds, 
        x_i=x_i,
        t=t,
        ndim=ndim,
        s_x=s_x,
        lhs_mat=lhs_mat,
        rhs_vec=rhs_vec,
        f_x=f_x,
        f_s_x=f_s_x,
        f=f
    )
    
    U, Σ, V = svd(lhs_mat)
    return V * ((transpose(U) * rhs_vec) ./ Σ), Σ
end


"""
    Used singular value decomposition to find the solution to`A Δx = g`, the linearised UPO problem. 

    # Arguments
    - `tds::TangentDynamicalSystem`: As in the `DynamicalSystems` library, tangent linear system of the model.
    - `x_i::Vector{AbstractFloat}`: Initial condition or initial guess for finding the UPO.
    - `t::AbstractFloat`: Initial guess of the period of the UPO.
    - `A_w::Matrix{AbstractFloat}`: Matrix defining the left hand side of the linear UPO problem `A Δx=g`.
    - `V::Adjoint{AbstractFloat, Matrix{AbstractFloat}}`: The value `V*` in `U Σ V* = A`.
    - `ell1::Vector{AbstractFloat}`: The value `Δx` in the linear problem `A Δx=g`.
    - `δ::AbstractFloat`: Defines the magnitude of purturbations.
    - `ndim::Int`: Dimension of the system.
    
    # Returns
    `Matrix{Float64}`: The 9 possible permutations of derivatives of quadratic terms of the UPO newton minimisation step.

"""
function tensor_correction(;
        tds::TangentDynamicalSystem, 
        x_i::Vector{T}, 
        t::T,
        A_w::Matrix{T},
        V::Adjoint{T, Matrix{T}}, 
        ell1::Vector{T}, 
        δ::T,
        ndim::Int
    ) where {T<:AbstractFloat}

    v_ell = V * ell1

    # Integration of linear system
    v_basis = Matrix{Float64}(undef, ndim, 9)

    v_basis[:, 1:3] .= V[1:end-1, end-1]
    v_basis[:, 1:1] .= integrate_tgls(tds, x_i .+ δ * V[1:ndim, end-1], t=t + δ * V[end, end-1], basis=v_basis[:, 1:1], return_x=false)
    v_basis[:, 2:2] .= integrate_tgls(tds, x_i .+ δ * V[1:ndim, end], t=t + δ * V[end, end], basis=v_basis[:, 2:2], return_x=false)
    v_basis[:, 3:3] .= integrate_tgls(tds, x_i .+ δ * v_ell[1:ndim], t=t + δ * v_ell[end], basis=v_basis[:, 3:3], return_x=false)

    v_basis[:, 4:6] .= V[1:end-1, end]
    v_basis[:, 4:4] .= integrate_tgls(tds, x_i .+ δ * V[1:ndim, end-1], t=t + δ * V[end, end-1], basis=v_basis[:, 4:4], return_x=false)
    v_basis[:, 5:5] .= integrate_tgls(tds, x_i .+ δ * V[1:ndim, end], t=t + δ * V[end, end], basis=v_basis[:, 5:5], return_x=false)
    v_basis[:, 6:6] .= integrate_tgls(tds, x_i .+ δ * v_ell[1:ndim], t=t + δ * v_ell[end], basis=v_basis[:, 6:6], return_x=false)

    v_basis[:, 7:9] .= v_ell[1:end-1]
    v_basis[:, 7:7] .= integrate_tgls(tds, x_i .+ δ * V[1:ndim, end-1], t=t + δ * V[end, end-1], basis=v_basis[:, 7:7], return_x=false)
    v_basis[:, 8:8] .= integrate_tgls(tds, x_i .+ δ * V[1:ndim, end], t=t + δ * V[end, end], basis=v_basis[:, 8:8], return_x=false)
    v_basis[:, 9:9] .= integrate_tgls(tds, x_i .+ δ * v_ell[1:ndim], t=t + δ * v_ell[end], basis=v_basis[:, 9:9], return_x=false)


    # Construction of the derivatives
    A_vN = A_w * V[:, end-1]
    A_vN1 = A_w * V[:, end]
    A_d_l = A_w * ell1

    v_basis[:, 1:1] .= 1 / δ .* (v_basis[:, 1:1] .- A_vN[1:end-1])
    v_basis[:, 2:2] .= 1 / δ .* (v_basis[:, 2:2] .- A_vN[1:end-1])
    v_basis[:, 3:3] .= 1 / δ .* (v_basis[:, 3:3] .- A_vN[1:end-1])

    v_basis[:, 4:4] .= 1 / δ .* (v_basis[:, 4:4] .- A_vN1[1:end-1])
    v_basis[:, 5:5] .= 1 / δ .* (v_basis[:, 5:5] .- A_vN1[1:end-1])
    v_basis[:, 6:6] .= 1 / δ .* (v_basis[:, 6:6] .- A_vN1[1:end-1])

    v_basis[:, 7:7] .= 1 / δ .* (v_basis[:, 7:7] .- A_d_l[1:end-1])
    v_basis[:, 8:8] .= 1 / δ .* (v_basis[:, 8:8] .- A_d_l[1:end-1])
    v_basis[:, 9:9] .= 1 / δ .* (v_basis[:, 9:9] .- A_d_l[1:end-1])

    return v_basis
end


"""
    Inplace function to describe a system of two equations of the quadratic terms of the UPO newton minimisation step. 

    # Arguments
    - `x::Vector{AbstractFloat}`: Initial condition or initial guess .
    - `_ell::Vector{AbstractFloat}`: The value `Δx` in the linear problem `A Δx=g`.
    - `coefs::Matrix{AbstractFloat}`: The 9 possible permutations of derivatives of quadratic terms of the UPO newton minimisation step.
    - `g_prime:Vector{AbstractFloat}`: Right hand side vector defined in the linear problem `U* g`.
    - `U::Matrix{AbstractFloat}`: The value `U` in `U Σ V* = A`.
    - `Σ::Vector{AbstractFloat}`: The value `Σ` in `U Σ V* = A`.
    
    # Returns
    `nothing`

"""
function tensor_sub_system!(
        x::Vector{T},
        _ell::Vector{T},
        coefs::Matrix{T},
        g_prime::Vector{T}, 
        U::Matrix{T}, 
        Σ::Vector{T}
        ) where {T<:AbstractFloat}

    nl_terms_1 = ((U[1:end-1, end-1] ⋅ coefs[:, 1]) * _ell[1] * _ell[1]
                + (U[1:end-1, end-1] ⋅ coefs[:, 5]) * _ell[2] * _ell[2]
                + (U[1:end-1, end-1] ⋅ (coefs[:, 4] .+ coefs[:, 2])) * _ell[1] * _ell[2]
                + (U[1:end-1, end-1] ⋅ (coefs[:, 3] .+ coefs[:, 7])) * _ell[1]
                + (U[1:end-1, end-1] ⋅ (coefs[:, 8] .+ coefs[:, 6])) * _ell[2]
                + (U[1:end-1, end-1] ⋅ coefs[:, 9]))

    nl_terms_2 = ((U[1:end-1, end] ⋅ coefs[:, 1]) * _ell[1] * _ell[1]
                + (U[1:end-1, end] ⋅ coefs[:, 5]) * _ell[2] * _ell[2]
                + (U[1:end-1, end] ⋅ (coefs[:, 4] .+ coefs[:, 2])) * _ell[1] * _ell[2]
                + (U[1:end-1, end] ⋅ (coefs[:, 3] .+ coefs[:, 7])) * _ell[1]
                + (U[1:end-1, end] ⋅ (coefs[:, 8] .+ coefs[:, 6])) * _ell[2]
                + (U[1:end-1, end] ⋅ coefs[:, 9]))
    
    x[1] = Σ[end-1] * _ell[1] + 1 / 2 * nl_terms_1 - g_prime[end-1]
    x[2] = Σ[end] * _ell[2] + 1 / 2 * nl_terms_2 - g_prime[end]
    return nothing

end


"""
    Finds the roots of the subsystem defined by `tensor_sub_system!`. 

    # Arguments
    - `subsystem::Function`: Describes a system of two equations of the quadratic terms of the UPO newton minimisation step. 
    - `max_its::Int`: Maximum number of attempts to find a minimum starting from a random initial condition.
    - `coefs::Matrix{AbstractFloat}`: The 9 possible permutations of derivatives of quadratic terms of the UPO newton minimisation step.
    - `g_prime:Vector{AbstractFloat}`: Right hand side vector defined in the linear problem `U* g`.
    - `U::Matrix{AbstractFloat}`: The value `U` in `U Σ V* = A`.
    - `Σ::Vector{AbstractFloat}`: The value `Σ` in `U Σ V* = A`.
    
    # Returns
    `Vector{AbstractFloat}`: Values that minimises `tensor_sub_system!`.

"""
function roots_of_tensor_correction(;
        subsystem::Function, 
        max_its::Int,
        coefs::Matrix{T},
        g_prime::Vector{T}, 
        U::Matrix{T}, 
        Σ::Vector{T}
    ) where {T<:AbstractFloat}

    found_roots = zeros(max_its, 2)
    
    func_closure = (err, x) -> subsystem(err, x, coefs, g_prime, U, Σ)

    iterations = 1
    total_loops = 0
    while (iterations < max_its+1) & (total_loops < max_its * 1e3)
        ic = randn(2) * 10
        r_test = nlsolve(func_closure, ic)
        if r_test.f_converged
            found_roots[iterations, :] .= r_test.zero
            iterations += 1
        end
        total_loops += 1
    end

    if iterations == 1
        output = [0., 0.]
    else
        error = Array{Float64}(undef, iterations-1)

        for i in 1:iterations-1
            x = [0., 0.]
            func_closure(x, found_roots[i, :])
            error[i] = norm(x)
        end
        output = found_roots[argmin(error), :]
    end
    return output 
end


"""
    Used singular value decomposition to find the solution to`A Δx = g`, the linearised UPO problem. 

    # Arguments
    - `tds::TangentDynamicalSystem`: As in the `DynamicalSystems` library, tangent linear system of the model.
    - `tds_one::TangentDynamicalSystem`: As in the `DynamicalSystems` library, tangent linear system of the model with only a single basis vector that is integrated.
    - `x_i::Vector{AbstractFloat}`: Initial condition or initial guess for finding the UPO.
    - `t::AbstractFloat`: Initial guess of the period of the UPO.
    - `ndim::Int`: Dimension of the system.
    - `s_x::Vector{AbstractFloat}`: Vector of the difference between the initial guess `x_i`, and the point in the phase space after integration time of `t` from the initial point `x_i`.
    - `lhs_mat::Matrix{AbstractFloat}`: Matrix `A` that describes the linearised UPO problem. This variable is over written so any matrix of the correct dimension can be passed.
    - `rhs_vec::Vector{AbstractFloat}`: Vector `g` that describes the linearised UPO problem. This variable is over written so any vector of the correct dimension can be passed.
    - `f_x::Array{AbstractFloat}`: The tendencies evaluated at the initial point `x_i`, this variable is over written so any vector of the correct dimension can be passed.
    - `f_s_x::Array{AbstractFloat}`: The tendencies evaluated at the point in phase space defined by the initial point `x_i` integrated for time `t`, this variable is over written so any vector of the correct dimension can be passed.
    - `f::Function`: Tendencies of the model `dx/dt=f(x)`.
    - `δ::AbstractFloat`: Defines the magnitude of purturbations.
    - `max_its::Int`: Maximum number of attempts to find a minimum starting from a random initial condition.
    
    # Returns
    `Vector{AbstractFloat}`: Correction step for the quadratic Newton step for finding UPOs.

"""
function tensor_problem(;
        tds::TangentDynamicalSystem,
        tds_one::TangentDynamicalSystem, 
        x_i::Vector{T},
        t::T,
        ndim::Int,
        s_x::Vector{T},
        lhs_mat::Matrix{T},
        rhs_vec::Vector{T},
        f_x::Array{T},
        f_s_x::Array{T},
        f::Function,
        delta::T,
        max_its::Int
    ) where {T<:AbstractFloat}
    
    # Calculate A (lhs matrix), and g (rhs vector)
    construct_matrix!(
        tds=tds, 
        x_i=x_i,
        t=t,
        ndim=ndim,
        s_x=s_x,
        lhs_mat=lhs_mat,
        rhs_vec=rhs_vec,
        f_x=f_x,
        f_s_x=f_s_x,
        f=f
    )
    
    U, Σ, V = svd(lhs_mat)
    g_prime = transpose(U) * rhs_vec

    # Solving the linear system
    ell1_prime = zeros(Float64, ndim + 1)
    # ell1_prime[1:end-2] .= transpose(U)[1:end-2, 1:end-2] * rhs_vec[1:end-2] ./ Σ[1:end-2]
    ell1_prime[1:end-2] .= Σ[1:end-2].^(-1) .* g_prime[1:end-2]

    # Calculating the 2-dim subsystem
    coefs = tensor_correction(
                tds=tds_one, 
                x_i=x_i,
                t=t, 
                A_w=lhs_mat, 
                V=V, 
                ell1=ell1_prime, 
                δ=delta,
                ndim=ndim)
    
    # Calculating the roots of the subsystem
    roots_of_correction = roots_of_tensor_correction(
                            subsystem=tensor_sub_system!, 
                            max_its=max_its,
                            coefs=coefs, 
                            g_prime=g_prime, 
                            U=U, 
                            Σ=Σ)

    ell1_prime[ndim:end] .= roots_of_correction
    return V * ell1_prime
end


"""
    Line search for finding the optimal step size in the Newton method.

    # Arguments
    - `ds::CoupledODEs`: As in the `DynamicalSystems` library, system of model equations. 
    - `x_i::Vector{AbstractFloat}`: Initial condition or initial guess for finding the UPO.
    - `s_x::Vector{AbstractFloat}`: The point in the phase space after integration time of `t` from the initial point `x_i`.
    - `t::AbstractFloat`: Initial guess of the period of the UPO.
    - `Δx::Vector{AbstractFloat}`: Vector of the correction step.
    - `ndim::Int`: Dimension of the system.
    - `del_x::Vector{AbstractFloat}`: Vector of the difference between the initial guess `x_i`, and the point in the phase space after integration time of `t` from the initial point `x_i`.
    - `max_steps::Int`: Maximum number of iterations of the algorithm.
    - `c::AbstractFloat`: Coefficient of magnitude of initial error.
    - `α::AbstractFloat`: Coefficient of the step size.
    - `ρ::AbstractFloat`: Reduction in magnitude of `α` in each iteration.
    - `max_period::AbstractFloat`: Bound on period of UPO.
    - `obj_func::Function`: Objective function to provide a metric for the error. Defualts to the Euclidian norm.
    - `obj_func_grad::Function`: Derivative of `obj_func`.

    # Returns
    `AbstractFloat`: Coefficient of the step size to apply to the error step.
    `Int`: The number of iterations used to determine this coefficient.
"""
function line_search(;
        ds::CoupledODEs, 
        x_i::Vector{T}, 
        s_x::Vector{T},
        t::T,
        Δx::Vector{T},
        ndim::Int,
        del_x::Vector{T},
        max_steps::Int,
        c::T,
        α::T,
        ρ::T,
        max_period::T,
        obj_func::Function,
        obj_func_grad::Function
    ) where {T<:AbstractFloat}

    integrate_traj!(ds, x_i, s_x, t)

    del_x[1:ndim] .= s_x - x_i
    del_x[end] = Δx[end]
    f_x_0 = obj_func(del_x)

    check = true
    k = 0
    r = c * dot(obj_func_grad(del_x), Δx)
    _α = copy(α)

    while check & (k < max_steps) & ((t + _α .* Δx[end]) < max_period)

        integrate_traj!(ds, x_i .+ _α .* Δx[1:end-1], s_x, max(t + _α .* Δx[end], 1e-5))

        del_x[1:end-1] .= s_x .- (x_i .+ _α .* Δx[1:end-1])
        del_x[end] = _α .* Δx[end]

        if (f_x_0 - obj_func(del_x)) >= _α * r
            check = false
        else
            k += 1
            _α = _α * ρ
        end
    end

    return _α, k
end


"""
    Hidden function that conducts one iteration of the Newton method for finding UPOs.
    
    # Argument
    - `tds::TangentDynamicalSystem`: As in the `DynamicalSystems` library, tangent linear system of the model.
    - `tds_one::TangentDynamicalSystem`: As in the `DynamicalSystems` library, tangent linear system of the model with only a single basis vector that is integrated.
    - `ds::CoupledODEs`: As in the `DynamicalSystems` library, the ODEs that describe the system.
    - `x_i::Vector{AbstractFloat}`: Initial condition or initial guess for finding the UPO.
    - `ndim::Int`: Dimension of the system.
    - `s_x::Vector{AbstractFloat}`: Vector of the difference between the initial guess `x_i`, and the point in the phase space after integration time of `t` from the initial point `x_i`.
    - `lhs_mat::Matrix{AbstractFloat}`: Matrix `A` that describes the linearised UPO problem. This variable is over written so any matrix of the correct dimension can be passed.
    - `rhs_vec::Vector{AbstractFloat}`: Vector `g` that describes the linearised UPO problem. This variable is over written so any vector of the correct dimension can be passed.
    - `f_x::Array{AbstractFloat}`: The tendencies evaluated at the initial point `x_i`, this variable is over written so any vector of the correct dimension can be passed.
    - `f_s_x::Array{AbstractFloat}`: The tendencies evaluated at the point in phase space defined by the initial point `x_i` integrated for time `t`, this variable is over written so any vector of the correct dimension can be passed.
    - `f::Function`: Tendencies of the model `dx/dt=f(x)`.
    - `damping::Bool`: Switch for the line search algorithm.
    - `damping_α::AbstractFloat`: Coefficient of the step size in the line search algorithm.
    - `damping_c::AbstractFloat`: Coefficient of magnitude of initial error in the line search algorithm.
    - `damping_ρ::AbstractFloat: Reduction in magnitude of `α` in each iteration in the line search algorithm.
    - `damping_max_steps::Int`: Maximum number of iterations of the line search algorithm.
    - `min_b::Union{Nothing, AbstractFloat}`: Minimum bound on each value of the UPO initial condition.
    - `max_b::Union{Nothing, AbstractFloat}`: Maximum bound on each value of the UPO initial condition
    - `max_period::AbstractFloat`: Initial guess of the period of the UPO.
    - `tensor_correction::Bool`: Switch to determine if the tensor correction method is used or not.
    - `tensor_auto::Bool`: Switch if the tensor method is automatically used, on condition of the magnitude of the singular values of the linear problem.
    - `tensor_trigger::AbstractFloat`: Threshold of the magnitude of the singular values of the linear problem, when a singular value is smaller than this threshold the tensor correction method is used.
    - `tensor_δ::AbstractFloat`: Defines the magnitude of purturbations in `roots_of_tensor_correction`.
    - `tensor_max_its::Int`: Maximum number of iterations of the tensor correction algorithm.

    # Returns
    - `AbstractFloat`: Resulting point after the correction step.
    - `AbstractFloat`: Error between the initial point and the corrected point.
    - `Int`: Number of line search steps used.
    - `String`: Type of algorithm used to make correction step.

"""
function _iterate_problem(;
        tds::TangentDynamicalSystem,
        tds_one::TangentDynamicalSystem,
        ds::CoupledODEs,
        x_i::Vector{T},
        ndim::Int,
        s_x::Vector{T},
        lhs_mat::Matrix{T},
        rhs_vec::Vector{T},
        f_x::Array{T},
        f_s_x::Array{T},
        f::Function,
        damping::Bool,
        damping_α::T,
        damping_c::T,
        damping_ρ::T,
        damping_max_steps::Int,
        min_b::Union{Nothing, T},
        max_b::Union{Nothing, T},
        max_period::T,
        tensor_correction::Bool,
        tensor_auto::Bool,
        tensor_trigger::T,
        tensor_δ::T,
        tensor_max_its::Int
    ) where {T<:AbstractFloat}
    
    if tensor_correction
        Δx = tensor_problem(
            tds=tds,
            tds_one=tds_one, 
            x_i=x_i[1:ndim],
            t=x_i[end],
            ndim=ndim,
            s_x=s_x,
            lhs_mat=lhs_mat,
            rhs_vec=rhs_vec,
            f_x=f_x,
            f_s_x=f_s_x,
            f=f,
            delta=tensor_δ,
            max_its=tensor_max_its
        )
        solution_method = "tensor"
    else

        Δx, Σ = svd_problem(
            tds=tds, 
            x_i=x_i[1:ndim], 
            t=x_i[end], 
            ndim=ndim, 
            s_x=s_x,
            lhs_mat=lhs_mat,
            rhs_vec=rhs_vec,
            f_x=f_x,
            f_s_x=f_s_x,
            f=f)
        
        solution_method = "linear"

        if tensor_auto
            if Σ[end] < tensor_trigger
                Δx = tensor_problem(
                    tds=tds,
                    tds_one=tds_one, 
                    x_i=x_i[1:ndim],
                    t=x_i[end],
                    ndim=ndim,
                    s_x=s_x,
                    lhs_mat=lhs_mat,
                    rhs_vec=rhs_vec,
                    f_x=f_x,
                    f_s_x=f_s_x,
                    f=f,
                    delta=tensor_δ,
                    max_its=tensor_max_its
                )
                solution_method = "tensor"
            end
        end
    end

    η = damping_α
    k = 0

    if damping & (norm(Δx) < norm(x_i[1:ndim]) * 1e2)
        η, k = line_search(
            ds=ds, 
            x_i=x_i[1:ndim], 
            s_x=s_x,
            t=x_i[end],
            Δx=Δx, 
            ndim=ndim,
            del_x=rhs_vec,
            max_steps=damping_max_steps,
            c=damping_c,
            α=damping_α,
            ρ=damping_ρ,
            max_period=max_period,
            obj_func=objective_function,
            obj_func_grad=grad_obj_function
        )
    end

    ic_n = x_i + Δx .* η

    # Bounds on the point
    if !isnothing(min_b)
        ic_n = max.(ic_n, min_b)
    end
    if !isnothing(max_b)
        ic_n = min.(ic_n, max_b)
    end

    # Bound on the minimum period
    ic_n[end] = max(1e-10, ic_n[end])

    Δ = norm(ic_n - x_i)

    return ic_n, Δ, k, solution_method
end


"""
    Hidden function that controls the Newton method UPO search.
    
    # Argument
    - `tds::TangentDynamicalSystem`: As in the `DynamicalSystems` library, tangent linear system of the model.
    - `ds::CoupledODEs`: As in the `DynamicalSystems` library, the ODEs that describe the system.
    - `ic::Vector{AbstractFloat}`: Initial condition or initial guess for finding the UPO.
    - `iterations::Int`: Maximum number of iterations of the correction steps.
    - `print_report::Bool`: Switch to determine whether a report is printed after algorithm completion.
    - `bounds::Union{Nothing, Dict}`: Bounds on the initial condition.
    - `min_norm::AbstractFloat`: Minimum size of the error, errors under this threshold will stop the algorithm.
    - `error_flex::AbstractFloat`: Multiple of the previous error that the magnitude next error can be. This aims to allow the algorithm to pass through local minima. 
    - `damping::Bool`: Switch for the line search algorithm.
    - `damping_α::AbstractFloat`: Coefficient of the step size in the line search algorithm.
    - `damping_c::AbstractFloat`: Coefficient of magnitude of initial error in the line search algorithm.
    - `damping_ρ::AbstractFloat: Reduction in magnitude of `α` in each iteration in the line search algorithm.
    - `damping_max_steps::Int`: Maximum number of iterations of the line search algorithm.
    - `max_period::AbstractFloat`: Initial guess of the period of the UPO.
    - `tensor_correction::Bool`: Switch to determine if the tensor correction method is used or not.
    - `tensor_auto::Bool`: Switch if the tensor method is automatically used, on condition of the magnitude of the singular values of the linear problem.
    - `tensor_trigger::AbstractFloat`: Threshold of the magnitude of the singular values of the linear problem, when a singular value is smaller than this threshold the tensor correction method is used.
    - `tensor_δ::AbstractFloat`: Defines the magnitude of purturbations in `roots_of_tensor_correction`.
    - `tensor_max_its::Int`: Maximum number of iterations of the tensor correction algorithm.

    # Returns
    - `Array{AbstractFloat}`: Concatination of the point in state space with the period of the UPO.

"""
function _find_upo_iterations_nm(;
        tds::TangentDynamicalSystem,
        ds::CoupledODEs,
        ic::Vector{T}, 
        iterations::Int, 
        print_report::Bool, 
        bounds::Union{Nothing, Dict}, 
        min_norm::T,
        error_flex::T,
        damping::Bool,
        damping_α::T,
        damping_c::T,
        damping_ρ::T,
        damping_max_steps::Int,
        max_period::T,
        tensor_correction::Bool,
        tensor_auto::Bool,
        tensor_trigger::T,
        tensor_δ::T,
        tensor_max_its::Int
    ) where {T<:AbstractFloat}

    ndim = dimension(ds)

    if length(ic)-1 != ndim
        println("IC needs to include period as well.")
    end

    delta_x = Array{Float64}(undef, iterations)

    # Upper and lower bounds for each variable
    min_arr, max_arr = nothing, nothing
    if !isnothing(bounds)
        if "min" in keys(bounds)
            min_b = bounds["min"]
            min_arr = ones(length(ic)) * -Inf
            min_arr = construct_bounds(min_arr, min_b)

        end
        if "max" in keys(bounds)
            max_b = bounds["max"]
            max_arr = ones(length(ic)) * Inf
            max_arr = construct_bounds(max_arr, max_b)
        end
    else
        max_arr = nothing
        min_arr = nothing
    end

    xs = Array{Float64}(undef, iterations+1, length(ic))

    xs[1, :] = ic
    report = Array{Any}(undef, iterations+1, 6)

    min_step = true
    i = 1
    prev_error = norm(F_x(ds, ic))
    next_error = prev_error

    # Create re-used variables

    # One vector TDS
    tds_one = _construct_tds_one_vec(ds, tds)
    
    f_rule = dynamic_rule(ds)
    s_x, lhs_mat, rhs_vec, f_x, f_s_x = _prealocate_variables(ndim)
    del_x = ic[end]

    while (i <= iterations) & min_step & (prev_error * error_flex >= next_error) & (del_x[end] < max_period)
        del_x, dx, k, method = _iterate_problem(
            tds=tds,
            tds_one=tds_one,
            ds=ds,
            x_i=copy(xs[i, :]),
            ndim=ndim,
            s_x=s_x,
            lhs_mat=lhs_mat,
            rhs_vec=rhs_vec,
            f_x=f_x,
            f_s_x=f_s_x,
            f=f_rule,
            damping=damping,
            damping_α=damping_α,
            damping_c=damping_c,
            damping_ρ=damping_ρ,
            damping_max_steps=damping_max_steps,
            min_b=min_arr,
            max_b=max_arr,
            max_period=max_period,
            tensor_correction=tensor_correction,
            tensor_auto=tensor_auto,
            tensor_trigger=tensor_trigger,
            tensor_δ=tensor_δ,
            tensor_max_its=tensor_max_its
        )

        xs[i+1, :] .= del_x
        delta_x[i] = dx
        
        prev_error = next_error
        next_error = norm(F_x(ds, del_x))
        
        report[i, :] = [i, del_x[end], dx, next_error, k, method]
        i += 1
        if dx < min_norm
            min_step = false
        end
    end

    if print_report
        end_row = max(i-2, 1)
        pretty_table(report[1:end_row, :], header=["Iteration", "Period", "Step Size", "Error", "# Line search", "Method"])
    end
    
    if (i - 1) < 1
        println("Iterations not successful. This can be because the error is larger in the next step, or the maximum period has been exceeded.")
    end
    return xs[max(i-1, 1), :]
end


"""
    Newton method UPO search.
    
    # Argument
    - `system::System`: System of equations the describe the model. 
    - `ic::Vector{AbstractFloat}`: Initial condition or initial guess for finding the UPO.
    - `iterations::Int`: Maximum number of iterations of the correction steps.
    - `print_report::Bool`: Switch to determine whether a report is printed after algorithm completion.
    - `bounds::Union{Nothing, Dict}`: Bounds on the initial condition.
    - `min_norm::AbstractFloat`: Minimum size of the error, errors under this threshold will stop the algorithm.
    - `error_flex::AbstractFloat`: Multiple of the previous error that the magnitude next error can be. This aims to allow the algorithm to pass through local minima. 
    - `damping::Bool`: Switch for the line search algorithm.
    - `damping_α::AbstractFloat`: Coefficient of the step size in the line search algorithm.
    - `damping_c::AbstractFloat`: Coefficient of magnitude of initial error in the line search algorithm.
    - `damping_ρ::AbstractFloat: Reduction in magnitude of `α` in each iteration in the line search algorithm.
    - `damping_max_steps::Int`: Maximum number of iterations of the line search algorithm.
    - `max_period::AbstractFloat`: Initial guess of the period of the UPO.
    - `tensor_correction::Bool`: Switch to determine if the tensor correction method is used or not.
    - `tensor_auto::Bool`: Switch if the tensor method is automatically used, on condition of the magnitude of the singular values of the linear problem.
    - `tensor_trigger::AbstractFloat`: Threshold of the magnitude of the singular values of the linear problem, when a singular value is smaller than this threshold the tensor correction method is used.
    - `tensor_δ::AbstractFloat`: Defines the magnitude of purturbations in `roots_of_tensor_correction`.
    - `tensor_max_its::Int`: Maximum number of iterations of the tensor correction algorithm.

    # Returns
    - `UPO_sol`: Concatination of the point in state space with the period of the UPO.

"""
function find_upo_nm(
        system::System,
        ic::Union{Vector{T}, UPO_sol{T}};
        iterations::Int=50, 
        print_report::Bool=false, 
        bounds::Union{Nothing, Dict}=nothing, 
        min_norm::T=1e-10,
        error_flex::T=1.,
        damping::Bool=true,
        damping_α::T=0.9,
        damping_c::T=1.,
        damping_ρ::T=0.9,
        damping_max_steps::Int=100,
        max_period::T=500.,
        tensor_correction::Bool=false,
        tensor_auto::Bool=false,
        tensor_trigger::T=1e-2,
        tensor_δ::T=1e-4,
        tensor_max_its::Int=100
    ) where {T<:AbstractFloat}

    # CoupledODEs and TangentLinear system
    ds, tds = generate_ds(system)

    if ic isa UPO_sol
        ndim = system.ndim
        _ic = Vector{Float64}(undef, ndim+1)
        _ic[1:ndim] = ic.ic
        _ic[end] = ic.period
    else
        _ic = ic
    end
    
    res = _find_upo_iterations_nm(
        tds=tds,
        ds=ds,
        ic=_ic, 
        iterations=iterations, 
        print_report=print_report, 
        bounds=bounds, 
        min_norm=min_norm,
        error_flex=error_flex,
        damping=damping,
        damping_α=damping_α,
        damping_c=damping_c,
        damping_ρ=damping_ρ,
        damping_max_steps=damping_max_steps,
        max_period=max_period,
        tensor_correction=tensor_correction,
        tensor_auto=tensor_auto,
        tensor_trigger=tensor_trigger,
        tensor_δ=tensor_δ,
        tensor_max_its=tensor_max_its
    )

    # return res
    prior_error = norm(F_x(ds, _ic))
    final_error = norm(F_x(ds, res))
    
    suc = (prior_error > final_error)
    Λ = floquet_multipliers(tds, res[1:end-1], res[end])
    return UPO_sol(
        success=suc, 
        ic=res[1:end-1], 
        period=res[end], 
        error=final_error, 
        system=SystemDef(system),
        floquet_multipliers=Λ
        ) 
end
