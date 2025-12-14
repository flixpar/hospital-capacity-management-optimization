using JuMP

"""
    transfers_subproblem(model, arrivals, L, params)

Add transfer decisions and occupancy recursion, returning the model and the
transfer-related objective expression.
"""
function transfers_subproblem(model, arrivals, L, params)
    N = model[:N]
    T = model[:T]
    Topt = model[:Topt]

    # decision variables
    @variable(model, transfers[1:N, 1:N, Topt] ≥ 0, integer=params.integer)

    # auxillary
    @expression(model, transfers_[i=1:N, j=1:N, t=1:T], (t ∈ Topt) ? transfers[i,j,t] : 0)

    # occupancy
    @expression(model, admissions[i=1:N, t=1:T], arrivals[i,t] + sum(transfers_[:,i,t]) - sum(transfers_[i,:,t]))
    @expression(model, discharges[i=1:N, t=1:T], dot(admissions[i,1:t], L[i,t:-1:1]))
    @expression(model, occupancy[i=1:N, t=1:T], sum(admissions[i,1:t]) - sum(discharges[i,1:t]))

    # transfer constraint
    @constraint(model, [i=1:N, t in Topt], sum(transfers[i,:,t]) ≤ arrivals[i,t])

    # transfer budgets
    add_transfer_budget!(model, params.budgets)

    obj = compute_transfer_objective(model, params)
    return model, obj
end

function compute_transfer_objective(model, params)
    objective = @expression(model, AffExpr(0))

    N = model[:N]
    T = length(model[:Topt])

    if any(params.costs .!= 0)
        transfers = model[:transfers].data
        add_to_expression!(objective, dot(params.costs, sum(transfers, dims=3)))
    end

    if params.transfer_smoothness > 0
        transfers = model[:transfers].data
        @variable(model, transfer_smoothness_penalties[i=1:N,j=i:N,1:(T-1)] ≥ 0)
        @constraint(model, [i=1:N, j=i:N, t=1:(T-1)], transfer_smoothness_penalties[i,j,t] >= transfers[i,j,t] - transfers[i,j,t+1])
        @constraint(model, [i=1:N, j=i:N, t=1:(T-1)], transfer_smoothness_penalties[i,j,t] >= transfers[i,j,t+1] - transfers[i,j,t])
        add_to_expression!(objective, params.transfer_smoothness * sum(transfer_smoothness_penalties))
    end

    if params.occupancy_smoothness > 0
        occupancy = model[:occupancy].data
        @variable(model, occupancy_smoothness_penalties[1:N,1:(T-1)] ≥ 0)
        @constraint(model, [i=1:N, t=1:(T-1)], occupancy_smoothness_penalties[i,t] >= occupancy[i,t] - occupancy[i,t+1])
        @constraint(model, [i=1:N, t=1:(T-1)], occupancy_smoothness_penalties[i,t] >= occupancy[i,t+1] - occupancy[i,t])
        add_to_expression!(objective, params.occupancy_smoothness * sum(occupancy_smoothness_penalties))
    end

    if params.admissions_smoothness > 0
        admissions = model[:admissions].data
        @variable(model, admissions_smoothness_penalties[1:N,1:(T-1)] ≥ 0)
        @constraint(model, [i=1:N, t=1:(T-1)], admissions_smoothness_penalties[i,t] >= admissions[i,t] - admissions[i,t+1])
        @constraint(model, [i=1:N, t=1:(T-1)], admissions_smoothness_penalties[i,t] >= admissions[i,t+1] - admissions[i,t])
        add_to_expression!(objective, params.admissions_smoothness * sum(admissions_smoothness_penalties))
    end

    return objective
end

function add_transfer_budget!(model, transferbudget)
    N = model[:N]
    Topt = model[:Topt]
    transfers = model[:transfers]

    notinf(x) = !isinf(x)
    if any(notinf.(transferbudget.perhospitalpair))
        @constraint(model, [i=1:N, j=1:N, t in Topt], transfers[i,j,t] ≤ transferbudget.perhospitalpair[i,j])
    end
    if any(notinf.(transferbudget.perhospital))
        @constraint(model, [i=1:N, t in Topt], sum(transfers[i,:,t]) ≤ transferbudget.perhospital[i])
    end
    if notinf(transferbudget.total)
        @constraint(model, [i=1:N, j=1:N, t in Topt], sum(transfers) ≤ transferbudget.total)
    end

    return model
end
