include("system.jl")
include("base.jl")

using DynamicalSystems
using Distances


begin
    Base.@kwdef mutable struct UPO_IC_kwargs{T<:AbstractFloat}
        min_threshold::T=1e-1
        max_time::T=1e2
        step_size::Int=1
        Δt::T=0.01
        error_threshold::T=0.5
        filter_gap::T=0.1
        Δperiod_mult::T=1e-5
        Δnorm_error::T=1e-1
    end
end


"""
    Filters a collection of `UPO_sol` to return a collection where the periods are all at least `gap` apart.

    # Arguments
    - `upo_guesses::Union{Vector{Any}, Vector{UPO_sol}}`: Collection of `UPO_sol` to filter.
    - `gap::Float64`: Gap in the periods of the collection of `UPO_sol`
    
    # Returns
    - `Array{UPO_sol}`
"""
function filter_ic_guesses(
    upo_guesses::Union{Vector{Any}, Vector{UPO_sol}};
    gap::AbstractFloat=0.5
)
    if length(upo_guesses) == 0
        output = []
    else
        periods = [g.period for g in upo_guesses]
        errors = [g.error for g in upo_guesses]

        min_period, max_period = minimum(periods), maximum(periods)
        steps = max(round(Int, (max_period - min_period) / gap), 2)

        bins = LinRange(min_period, max_period, steps)
        
        output = []
        for (i, j) in zip(bins[1:end-1], bins[2:end])
            mask = (periods .>= i) .& (periods .<= j)
            es = errors[mask]
            if length(es) > 0
                ix_min = argmin(es)
                push!(output, upo_guesses[mask][ix_min])
            end
        end
    end
    return output
end


"""
    Checks if two UPOs describe the same UPO but their periods differ by a multiple.

    # Arguments
    - `upo1::UPO_sol`: UPO number 1
    - `upo2::UPO_sol`: UPO to compare with the first.
    - `Δperiod_mult::Float64=1e-5`: Error permitted in the period up to a multiple.
    - `Δnorm_error::Float64=1e-1`: Maximum magnitude between points that is permitted.
    
    # Returns
    - `Bool`, `true` if the UPOs are the same up to a multiple in the period.
"""
function unique_multiple_period(
    upo1::UPO_sol,
    upo2::UPO_sol;
    Δperiod_mult::T=1e-5,
    Δnorm_error::T=1e-1
)::Bool where {T<:AbstractFloat}

    # Check which has largest period
    if upo2.period > upo1.period
        u1, u2 = upo1, upo2
    else
        u1, u2 = upo2, upo1
    end

    # Check for periodicity
    period_mult = u2.period / u1.period
    if (round(Int, period_mult) - period_mult < Δperiod_mult)
        traj1, _ = trajectory(ls, u1.period, u1.ic, Δt=0.01)
        traj2, _ = trajectory(ls, u2.period, u2.ic, Δt=0.01)

        norms = zeros(length(traj1))
        for i in 1:length(traj1)
            norms[i] = norm(traj2[1] .- traj1[i])
        end
        
        ix = argmin(norms)
        traj_new, time_new = trajectory(ls, u2.period, traj1[ix], Δt=0.01)
        
        norms = zeros(length(traj2))
        for i in 1:length(traj2)
            norms[i] = norm(traj2[i] .- traj_new[i])
        end
        return maximum(norms) > Δnorm_error
        
    else
        return true
    end
end


"""
    Filters a collection of `UPO_sol` to return a collection where the periods are all at least `gap` apart.

    # Arguments
    - `upo_guesses::Union{Vector{Any}, Vector{UPO_sol}}`: Collection of `UPO_sol` to filter.
    - `gap::Float64`: Gap in the periods of the collection of `UPO_sol`
    
    # Returns
    - `UPO_sol`
"""
function upo_similarity_filter_exact(
    ls1::System,
    upo1::UPO_sol,
    upo2::UPO_sol;
    Δnorm_error::T=0.005,
    Δnorm_error_mean::T=0.001,
    ΔT::T=0.001,
    ls2::Union{System, Nothing}=nothing
    )::Bool where {T<:AbstractFloat}
    # Function checks two UPOs to see if they are likely to result in the same UPO
    output = false

    period_error = abs(upo1.period - upo2.period)

    ls1 = generate_ode(ls1)

    if isnothing(ls2)
        ls2 = ls1
    else
        ls2 = generate_ode(ls2)
    end

    if abs(upo1.period - upo2.period) < ΔT
        traj1, _ = trajectory(ls1, upo1.period, upo1.ic, Δt=0.01)
        traj2, _ = trajectory(ls2, upo2.period, upo2.ic, Δt=0.01 * upo2.period / upo1.period)

        traj1 = Matrix(traj1)
        traj2 = Matrix(traj2)

        vec = Vector{Float64}(undef, length(traj1))
        start_2_ix = argmin(norm.(eachcol(transpose(traj1) .- traj2[1, :])))
        
        del = Vector{Float64}(undef, size(traj1)[1])

        if abs(size(traj1)[1] - size(traj2)[1]) < 1
            for i in 1:min(size(traj1)[1], size(traj2)[1])-1
                t2_ix = ((i+start_2_ix-1) % size(traj1)[1]) + 1
                del[i] = norm(traj1[i, :] .- traj2[t2_ix, :])
            end
        else
            del[1] = 1e10
            println("Error in array lengths")
        end
        
        max_norm = maximum(del)
        # println(max_norm)
        # println(mean(del))

        if (max_norm < Δnorm_error) | (mean(del) < Δnorm_error_mean)
            output = true
        end
    end
    return output
end


"""
    Removes UPOs in an array that are the same up to a multiple of the period.

    # Arguments
    - `upo_guesses::Union{Vector{Any}, Vector{UPO_sol}}`: Collection of `UPO_sol` to filter.
    - `Δperiod_mult::Float64=1e-5`: Error permitted in the period up to a multiple.
    - `Δnorm_error::Float64=1e-1`: Maximum magnitude between points that is permitted.
    
    # Returns
    - `Array{UPO_sol}`
"""
function remove_period_doublings(
        upo_guesses::Union{Vector{Any}, Vector{UPO_sol}}; 
        Δperiod_mult::T=1e-5,
        Δnorm_error::T=1e-1
    ) where {T<:AbstractFloat}
    upo_guesses_copy = copy(upo_guesses)
    upo_outputs = [upo_guesses[1]]
    deleteat!(upo_guesses_copy, 1)
    
    while length(upo_guesses_copy) > 0
        u2 = upo_guesses_copy[1]

        keep_upo = true
        i = 1
        while keep_upo & (i <= length(upo_outputs))
            keep_upo = unique_multiple_period(
                                upo_outputs[i], 
                                u2, 
                                Δperiod_mult=Δperiod_mult, 
                                Δnorm_error=Δnorm_error
                            )
            i += 1
        end

        if keep_upo
            push!(upo_outputs, u2)
        end
        deleteat!(upo_guesses_copy, 1)
    end

    return upo_outputs
end  


"""
    Finds near misses of a periodic orbit by following a trajectory and finding points where the trajectory passes close to the point in the future.

    # Arguments
    - `system::System`: System of equations the describe the model.
    - `max_period::Float64=2.`: Maximum period the algorithm looks for.
    - `ic::Union{Nothing, Vector{Float64}}=nothing`: Initial condition of the trajectory. If `nothing` a random initial condition is used.
    - `min_period::Float64=1.`: Minimum period the algorithm looks for. Minimum periods should not be set too small as the method will likely return points that are short lines in phase space rather than orbits.
    - `filter_guesses::Bool=true`: Filters the UPOs found so only one is returned of a given period. This uses `filter_ic_guesses`.
    - `filter_periods::Bool=false`: Removes periods that are the same up to mulitples in the period. This uses `remove_period_doublings`
    - `ic_kwargs::UPO_IC_kwargs=UPO_IC_kwargs()
    
    # Returns
    - `Array{UPO_sol}`
"""
function guess_ic(
        system::System;
        max_period::T=2.,
        ic::Union{Nothing, Vector{T}}=nothing,
        min_period::T=1.,
        filter_guesses::Bool=true,
        filter_periods::Bool=false,
        ic_kwargs::UPO_IC_kwargs=UPO_IC_kwargs()
    ) where {T<:AbstractFloat}

    bools = Vector{Bool}(undef, 0)
    period_save = Vector{Float64}(undef, 0)
    ics = Matrix{Float64}(undef, 0, system.ndim)
    Δ = Vector{Float64}(undef, 0)
    
    if isnothing(ic)
        ic = rand(system.ndim)
    end

    ls = generate_ode(system)

    traj, time = trajectory(ls, ic_kwargs.max_time, ic, Δt=ic_kwargs.Δt)
    max_steps = round(Int, max_period / ic_kwargs.Δt)
    min_steps = round(Int, min_period / ic_kwargs.Δt)

    end_step = length(traj) - max_steps
    # start_j = 
    
    for i in 1:ic_kwargs.step_size:end_step
        for j in i+min_steps:ic_kwargs.step_size:i + max_steps
            dist = norm(traj[i] - traj[j])
            if dist < ic_kwargs.error_threshold
                period = time[j] - time[i]
                if period > min_period
                    push!(bools, true)
                    push!(period_save, period)
                    ics = [ics; traj[i]']
                    push!(Δ, dist)
                end
            end
        end
    end

    output = []
    for i in 1:length(bools)
        if bools[i]
            _temp_sol = UPO_sol(
                success=bools[i],
                ic=ics[i, :],
                period=period_save[i],
                error=Δ[i],
                system=SystemDef(system)
            )

            push!(output, _temp_sol)
        end
    end

    if filter_guesses
        output = filter_ic_guesses(output, gap=ic_kwargs.filter_gap)
    end
    if filter_periods
        output = remove_period_doublings(
                            output, 
                            Δperiod_mult=ic_kwargs.Δperiod_mult,
                            Δnorm_error=ic_kwargs.Δnorm_error
        )
    end
    return output
end

