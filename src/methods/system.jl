include("../systems/toy_systems.jl")
include("../systems/MAOLAM2x2.jl")

struct System{T<:AbstractFloat}
    name::String
    f::Function
    jac::Function
    ndim::Int
    p::Vector{T}
    p_names::Vector{String}
end

Base.copy(s::System) = System(
                        s.name, 
                        s.f, 
                        s.jac, 
                        s.ndim, 
                        copy(s.p), 
                        s.p_names
                        )

struct SystemDef{T<:AbstractFloat}
    name::String
    ndim::Int
    p::Vector{T}
    p_names::Vector{String}
end

SystemDef(sys::System) = SystemDef(
    sys.name,
    sys.ndim,
    sys.p,
    sys.p_names
)

function system_to_dict(sys::SystemDef)
    return Dict(
        "name" => sys.name,
        "ndim" => sys.ndim,
        "p" => sys.p,
        "p_names" => sys.p_names
    )
end

function dict_to_system(d::Dict{String, Any})::SystemDef
    return SystemDef(
        d["name"],
        d["ndim"],
        map(x -> Float64(x), d["p"]),
        map(x -> String(x), d["p_names"])
    )
end

begin
    lorenz = System(
        "Lorenz", 
        lorenz!, 
        lorenz_j!, 
        3, 
        [10., 28., 8/3],
        ["sig", "r", "b"]
        )
end

begin
    maolam_2x2 = System(
        "MALOAM-2x2-heat_transfer", 
        MAOLAM2x2_tendencies!, 
        MAOLAM2x2_jac!, 
        30, 
        [300., 0.4 * 300., 0.085, 0.02, 0.76],
        ["Cg1", "Ca1", "k_d", "k_p", "eps"]
        )
end


function _p_name_to_p(sys::System, p_names::Union{String, Vector{String}})
    vals = []

    p_n = sys.p_names

    if typeof(p_names) == Vector{String}
        for p in p_names
            loc = findall(p_n -> p_n==p, p_n)
            if length(loc) > 1
                println("More than one parameter with the same name!")
            end
            push!(vals, loc[1])
        end
    else    
        loc = findall(p_n->p_n==p_names, p_n)
        if length(loc) > 1
            println("More than one parameter with the same name!")
        end
        push!(vals, loc[1])
    end
    return vals[1]
end