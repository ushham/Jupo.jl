include("system.jl")
include("base.jl")
include("newton_method.jl")


@enum UPOMethod begin
    initial_condition
    newton_method
    stabalising_transforms
    chaos_control
end


function upo_finder(
    ds::TangentDynamicalSystem, 
    ic::Vector{Float64};
    method::UPOMethod=newton_method,
    iterations=50, 
    print_report=false, 
    bounds=nothing, 
    min_norm=1e-10,
    tensor=Tensor_Options(),
    damping=Damping_Options(),
    ic_conditions=UPO_IC_kwargs()
    )
    
    if method == initial_condition
        output = guess_ic(ds, )

    elseif method == newton_method
        output = find_upo_nm(
            ds::TangentDynamicalSystem, 
            ic::Vector{Float64};
            iterations=50, 
            damping=Damping_Options(), 
            print_report=false, 
            bounds=nothing, 
            tensor=Tensor_Options(), 
            min_norm=1e-10
        )
    end

    return output

end