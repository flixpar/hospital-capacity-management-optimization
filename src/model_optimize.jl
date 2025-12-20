using JuMP
using Gurobi

include("los.jl")
include("model_capacity.jl")
include("model_transfers.jl")
include("validation.jl")

"""
    unpack_decisions(model)

Extract decision arrays from an optimized JuMP model.
"""
function unpack_decisions(model)
    N = model[:N]
    Topt = model[:Topt]
    B = model[:B]

    capacity = Array(value.(model[:capacity]))
    admissions = value.(model[:admissions])[:, Topt]
    occupancy = value.(model[:occupancy])[:, Topt]

    transfers = if haskey(model, :transfers)
        Array(value.(model[:transfers]))
    else
        zeros(Int, N, N, length(Topt))
    end

    capacity_unit_allocated_ = value.(model[:capacity_unit_usable])
    capacity_unit_allocated = [[capacity_unit_allocated_[i,t,b] for b in 1:B[i]] for i in 1:N, t in Topt]

    shortage = if haskey(model, :shortage)
        value.(model[:shortage])
    else
        nothing
    end

    return (; transfers, capacity, capacity_unit_allocated, admissions, occupancy, shortage)
end

"""
    set_model_params!(model, params)

Stash metadata on the model for later convenience.
"""
function set_model_params!(model, params)
    for (k, v) in pairs(params)
        model[k] = v
    end
    return model
end

"""
    optimize_decisions(arrivals, capacity, los, Topt, capacity_params, transfer_params, solver_params; nonsurge_occupancy=nothing, total_capacity=nothing)

Build and solve the joint capacity/transfer optimization model.

Optional kwargs for non-surge patient shortage penalty:
- `nonsurge_occupancy`: Non-surge patient census matrix [N, T]
- `total_capacity`: Total staffed capacity per hospital [N]
"""
function optimize_decisions(
    arrivals::Array{<:Real,2},
    capacity::Array,
    los::Array{<:Distribution,1},
    Topt::Array{Int,1},
    capacity_params,
    transfer_params,
    solver_params;
    nonsurge_occupancy=nothing,
    total_capacity=nothing,
)
    validate_optimize_inputs(
        arrivals,
        capacity,
        los,
        Topt,
        capacity_params,
        transfer_params,
        solver_params;
        nonsurge_occupancy=nonsurge_occupancy,
        total_capacity=total_capacity,
    )

    N, T = size(arrivals)
    B = [length(capacity[i]) for i in 1:N]

    L = discretize_los(los, N, T)

    model = Model(Gurobi.Optimizer)
    if !solver_params.verbose
        set_silent(model)
    end
    if solver_params.timelimit > 0
        set_time_limit_sec(model, solver_params.timelimit)
    end
    set_model_params!(model, (; N, T, B, Topt))

    objective = @expression(model, AffExpr(0))

    if transfer_params.optimize
        model, transfers_objective = transfers_subproblem(model, arrivals, L, transfer_params)
        objective += transfers_objective
    end
    if capacity_params.optimize
        model, capacity_objective = capacity_subproblem(model, arrivals, capacity, L, capacity_params;
                                                        nonsurge_occupancy=nonsurge_occupancy,
                                                        total_capacity=total_capacity)
        objective += capacity_objective
    end

    @objective(model, Min, objective)
    optimize!(model)

    return model
end
