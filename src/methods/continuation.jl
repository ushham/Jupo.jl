include("system.jl")
include("base.jl")
include("newton_method.jl")

using DynamicalSystems
using LinearAlgebra
using PrettyTables


"""
    Saves 
"""
function continuation_save(upo::UPO_sol; folder::Union{String, Nothing})
    name = upo.system.name * "_" * string(round(upo.period, sigdigits=5))
    for (n, v) in zip(upo.system.p_names, upo.system.p)
        name = name * "_" * n * string(v)
    end
    save_upo(upo, file_name=name, folder=folder)
end


function _cnt_upo(
    system1::System,
    start_upo::UPO_sol, 
    par::String, 
    p_bound::Function,
    step_size::T,
    max_its::Int,
    p_ix::Int,
    Δnorm_error::T,
    Δnorm_error_mean::T,
    ΔT_error::T,
    upo_error::T,
    error_flex::T,
    folder::Union{String, Nothing}
    ) where {T<:AbstractFloat}

    
    its = 0
    v = start_upo.system.p[p_ix]

    _step = step_size
    _upo = start_upo

    system1.p[p_ix] = v
    system2 = copy(system1)

    while p_bound(v) & (its < max_its)
        v_try = v + _step / (2^its)
        
        println("Running for parameter: " * par * "=" * string(v_try))
        
        system2.p[p_ix] = v_try

        upo_try = find_upo_nm(
                system2, 
                _upo, 
                iterations=100,
                damping=true,
                damping_α=1.,
                damping_max_steps=100,
                tensor_auto=true,
                tensor_trigger=1e-4,
                tensor_max_its=100,
                max_period= ΔT_error * 10,
                error_flex=error_flex,
                min_norm=1e-10
            )
        
        its += 1
        if upo_similarity_filter_exact(
            system1, 
            _upo, 
            upo_try, 
            Δnorm_error=Δnorm_error, 
            Δnorm_error_mean=Δnorm_error_mean,  
            ΔT=ΔT_error, 
            ls2=system2
            )
            if upo_try.error < upo_error
                continuation_save(upo_try, folder=folder)
                _upo = upo_try
                system1.p[p_ix] = system2.p[p_ix]

                its = 0
                v = v_try
            end
        end

    end

end


function continue_upo(
    system1::System,
    start_upo::UPO_sol, 
    par::String, 
    p_bounds::Tuple{T, T}; 
    step_size::T=0.1,
    max_its::Int=10,
    Δnorm_error::T=0.005,
    Δnorm_error_mean::T=0.001,
    ΔT_error::T=2.,
    upo_error::T=1e-6,
    error_flex::T=2.,
    folder::Union{String, Nothing}=nothing
    ) where {T<:AbstractFloat}

    p_ix = _p_name_to_p(system1, par)

    # Increasing the variable
    println("~~~~~~ ↑↑↑  Continuing towards upper bound  ↑↑↑ ~~~~~~")
    bound_upper(v) = (v < p_bounds[2])
    _cnt_upo(
        system1,
        start_upo, 
        par,
        bound_upper,
        step_size,
        max_its,
        p_ix,
        Δnorm_error,
        Δnorm_error_mean,
        ΔT_error,
        upo_error,
        error_flex,
        folder
        )

    # Decreasing the variable
    println("~~~~~~ ↓↓↓  Continuing towards lower bound  ↓↓↓ ~~~~~~")
    bound_lower(v) = (v > p_bounds[1])
    _cnt_upo(
        system1,
        start_upo, 
        par,
        bound_lower,
        -step_size,
        max_its,
        p_ix,
        Δnorm_error,
        Δnorm_error_mean,
        ΔT_error,
        upo_error,
        error_flex,
        folder
        )
    
end


