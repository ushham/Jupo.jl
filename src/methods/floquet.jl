
"""
    Calculates the Floquet Multipliers of an periodic orbit.

    # Arguments
    - `tds::TangentDynamicalSystem`: Tangent linear system, as defined in the `DynamicalSystems` library.
    - `ic::Vector{AbstractFloat}`: Initial condition in phase space of the periodic orbit
    - `period::AbstractFloat`: Period of the periodic orbit.

    # Returns
    - `Union{Vector{Complex}, Vector{AbstractFloat}}`: Floquet multipliers of the periodic orbit.
"""
function floquet_multipliers(
    tds::TangentDynamicalSystem,
    ic::Vector{T},
    period::T
    )::Union{Vector{T}, Vector{Complex{T}}} where {T<:AbstractFloat}
    reinit!(tds, ic, Q0=Matrix(1.0I, length(ic), length(ic)))
    step!(tds, period, true)
    M = current_state(tds.ds)[:, 2:end]
    Λ = eigvals(M)
    return reverse(sort(Λ, by=abs))
end


"""
    Calculates the Floquet Exponents of an periodic orbit, given the Floquet Multipliers.

    # Arguments
    - `Λ::Union{Vector{Complex{AbstractFloat}}}, UPO_sol}`: Floquet multipliers of a periodic orbit.
    - `T_p::AbstractFloat`: Period of the periodic orbit.

    # Returns
    - `Floquet_exp`: `struct` containing the real and complex components o the Floquet exponents.
"""
function floquet_exponents(
    Λ::Union{Vector{Complex{T}}, UPO_sol{T}};
    T_p::Union{T, Nothing}=nothing
    ) where {T<:AbstractFloat}
    if not(isnothing(T_p))
        period = T_p
    end
    if Λ isa UPO_sol
        period = Λ.T
        Λ = Λ.floquet_multipiers
    end
    exps = log.(Λ) .* (1 / period)
    return Floquet_exp(real.(exps), imag.(exps))
end
