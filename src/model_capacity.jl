using JuMP

include("config.jl")

"""
    capacity_subproblem(model, arrivals, capacity_data, L, params; nonsurge_occupancy=nothing, total_capacity=nothing)

Add capacity allocation decisions and constraints to the model and return the
augmented model plus the capacity-related objective expression.

Optional kwargs for non-surge patient shortage penalty:
- `nonsurge_occupancy`: Non-surge patient census matrix [N, T]
- `total_capacity`: Total staffed capacity per hospital [N]
"""
function capacity_subproblem(model, arrivals, capacity_data, L, params; nonsurge_occupancy=nothing, total_capacity=nothing)
    N = model[:N]
    T = model[:T]
    Topt = model[:Topt]
    B = model[:B]

    # occupancy (reuse if transfers already added it)
    if !haskey(model, :occupancy)
        @expression(model, discharges[i=1:N, t=1:T], dot(arrivals[i,1:t], L[i,t:-1:1]))
        @expression(model, occupancy[i=1:N, t=1:T], sum(arrivals[i,1:t]) - sum(discharges[i,1:t]))
        model[:admissions] = arrivals
    else
        occupancy = model[:occupancy]
    end

    # decision variables
    if params.enforce_leadtime
        @variable(model, capacity_unit_allocated[i=1:N, Topt, 1:B[i]], Bin)
        @variable(model, capacity_unit_usable[i=1:N, Topt, 1:B[i]], Bin)

        st = [[capacity_data[i][b].time_setup for b in 1:B[i]] for i in 1:N]
        bt = [[capacity_data[i][b].time_breakdown for b in 1:B[i]] for i in 1:N]
        beds = [[capacity_data[i][b].beds for b in 1:B[i]] for i in 1:N]
        tmin, tmax = extrema(Topt)

        @constraint(model, [i=1:N, t in Topt, b=1:B[i], z=0:st[i][b]], capacity_unit_usable[i,t,b] ≤ capacity_unit_allocated[i,max(tmin,t-z),b])
        @constraint(model, [i=1:N, t in Topt, b=1:B[i], z=0:bt[i][b]], capacity_unit_usable[i,t,b] ≤ capacity_unit_allocated[i,min(tmax,t+z),b])

        @expression(model, capacity_usable[i=1:N, t in Topt], dot(capacity_unit_usable[i,t,:], beds[i]))
    else
        @variable(model, capacity_unit_usable[i=1:N, Topt, 1:B[i]], Bin)
        @expression(model, capacity_unit_allocated[i=1:N, t in Topt, b=1:B[i]], capacity_unit_usable[i,t,b])
    end

    # capacity
    @expression(model, capacity[i=1:N, t=Topt], sum(capacity_unit_usable[i,t,b] * capacity_data[i][b].beds for b in 1:B[i]))
    @constraint(model, [i=1:N, t in Topt], capacity[i,t] ≥ occupancy[i,t])

    # maximum occupancy rate
    if params.max_occupancy < 1
        @constraint(model, [i=1:N, t in Topt], params.max_occupancy * capacity[i,t] ≥ occupancy[i,t])
    end

    # order constraint
    if params.ordered
        @constraint(model, [i=1:N, t in Topt, b=2:B[i]], capacity_unit_usable[i,t,b] ≤ capacity_unit_usable[i,t,b-1])
    end

    if params.baseline_always
        capacity_levels = [[b.capacity_level for b in h] for h in capacity_data]
        Bbaseline = [findall(==(1), h) for h in capacity_levels]
        @constraint(model, [i=1:N, t in Topt, b in Bbaseline[i]], capacity_unit_usable[i,t,b] == 1)
    end

    obj = compute_capacity_objective(model, capacity_data, params; nonsurge_occupancy=nonsurge_occupancy, total_capacity=total_capacity)
    return model, obj
end

function compute_capacity_objective(model, capacity_costs, params; nonsurge_occupancy=nothing, total_capacity=nothing)
    objective = @expression(model, AffExpr(0))

    N = model[:N]
    Topt = model[:Topt]
    B = model[:B]

    capacity = model[:capacity].data
    capacity_unit_allocated = model[:capacity_unit_allocated]
    occupancy = model[:occupancy]

    if params.costs_unitday
        for i in 1:N, b in 1:B[i]
            if capacity_costs[i][b].cost_per_unitday > 0
                add_to_expression!(objective, capacity_costs[i][b].cost_per_unitday * sum(capacity_unit_allocated[i,:,b]))
            end
        end
    end

    if params.costs_bedday
        error("costs_bedday not implemented")
    end

    # capacity setup
    if params.costs_setup
        @variable(model, capacity_unit_setup[i=1:N, Topt, 1:B[i]], Bin)
        @constraint(model, [i=1:N, t in Topt[2:end], b=1:B[i]], capacity_unit_setup[i,t,b] ≥ capacity_unit_allocated[i,t,b] - capacity_unit_allocated[i,t-1,b])
        @constraint(model, [i=1:N, t in Topt[1:1], b=1:B[i]], capacity_unit_setup[i,t,b] == 1)

        for i in 1:N, b in 1:B[i]
            add_to_expression!(objective, capacity_costs[i][b].cost_setup * sum(capacity_unit_setup[i,:,b]))
        end
    end

    # capacity breakdown
    if params.costs_breakdown
        @variable(model, capacity_unit_breakdown[i=1:N, Topt, 1:B[i]], Bin)
        @constraint(model, [i=1:N, t in Topt[2:end], b=1:B[i]], capacity_unit_breakdown[i,t,b] ≥ capacity_unit_allocated[i,t,b] - capacity_unit_allocated[i,t+1,b])
        @constraint(model, [i=1:N, t in Topt[1:1], b=1:B[i]], capacity_unit_breakdown[i,t,b] == 0)

        for i in 1:N, b in 1:B[i]
            add_to_expression!(objective, capacity_costs[i][b].cost_breakdown * sum(capacity_unit_breakdown[i,:,b]))
        end
    end

    # capacity converted (setup/breakdown)
    if params.costs_convert
        @variable(model, capacity_unit_converted[i=1:N, Topt, 1:B[i]], Bin)
        @constraint(model, [i=1:N, t in Topt[2:end], b=1:B[i]], capacity_unit_converted[i,t,b] ≥ capacity_unit_allocated[i,t,b] - capacity_unit_allocated[i,t-1,b])
        @constraint(model, [i=1:N, t in Topt[2:end], b=1:B[i]], capacity_unit_converted[i,t,b] ≥ capacity_unit_allocated[i,t-1,b] - capacity_unit_allocated[i,t,b])
        @constraint(model, [i=1:N, t in Topt[1:1], b=1:B[i]], capacity_unit_converted[i,t,b] == 1)

        for i in 1:N, b in 1:B[i]
            add_to_expression!(objective, capacity_costs[i][b].cost_convert * sum(capacity_unit_converted[i,:,b]))
        end
    end

    if params.capacity_smoothness > 0
        T = length(model[:Topt])
        @variable(model, capacity_smoothness_penalties[1:N,1:(T-1)] ≥ 0)
        @constraint(model, [i=1:N, t=1:(T-1)], capacity_smoothness_penalties[i,t] >= capacity[i,t] - capacity[i,t+1])
        @constraint(model, [i=1:N, t=1:(T-1)], capacity_smoothness_penalties[i,t] >= capacity[i,t+1] - capacity[i,t])
        add_to_expression!(objective, params.capacity_smoothness * sum(capacity_smoothness_penalties))
    end

    # Non-surge patient shortage penalty (optional)
    # Penalizes when total census (surge + non-surge) exceeds total staffed capacity
    if params.shortage_penalty > 0 && !isnothing(nonsurge_occupancy) && !isnothing(total_capacity)
        # shortage[i,t] = max(0, surge_occupancy[i,t] + nonsurge_occupancy[i,t] - total_capacity[i])
        @variable(model, shortage[i=1:N, t in Topt] >= 0)
        @constraint(model, [i=1:N, t in Topt],
            shortage[i,t] >= occupancy[i,t] + nonsurge_occupancy[i,t] - total_capacity[i])

        add_to_expression!(objective, params.shortage_penalty * sum(shortage))
    end

    return objective
end
