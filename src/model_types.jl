struct TransferBudgets
    perhospitalday::Array{Real,2}
    perhospitalpair::Array{Real,2}
    perday::Array{Real,1}
    perhospital::Array{Real,1}
    total::Real
end

struct TransferParams
    optimize::Bool
    integer::Bool
    costs::Array{Real,2}
    budgets::TransferBudgets
    transfer_smoothness::Real
    occupancy_smoothness::Real
    admissions_smoothness::Real
end

struct CapacityParams
    optimize::Bool
    costs_setup::Bool
    costs_breakdown::Bool
    costs_convert::Bool
    costs_bedday::Bool
    costs_unitday::Bool
    ordered::Bool
    baseline_always::Bool
    enforce_leadtime::Bool
    capacity_smoothness::Real
    max_occupancy::Real
end

struct SolverParams
    timelimit::Real
    verbose::Bool
end

function TransferBudgets(N, T; perhospitalday=nothing, perhospitalpair=nothing, perday=nothing, perhospital=nothing, total=nothing)
    perhospitalday = isnothing(perhospitalday) ? fill(Inf, (N, T)) : (ndims(perhospitalday) == 1) ? repeat(perhospitalday, 1, T) : perhospitalday
    perhospitalpair = isnothing(perhospitalpair) ? fill(Inf, (N, N)) : perhospitalpair
    perday = isnothing(perday) ? fill(Inf, T) : (perday isa Array) ? perday : fill(perday, T)
    perhospital = isnothing(perhospital) ? fill(Inf, N) : perhospital
    total = isnothing(total) ? Inf : total
    return TransferBudgets(perhospitalday, perhospitalpair, perday, perhospital, total)
end

function TransferParams(N, T; optimize=true, integer=true, costs=nothing, budgets=nothing, transfer_smoothness=0, occupancy_smoothness=0, admissions_smoothness=0)
    costs = isnothing(costs) ? fill(0, (N, N)) : (costs isa Array) ? costs : fill(costs, (N, N))
    budgets = isnothing(budgets) ? TransferBudgets(N, T) : budgets
    return TransferParams(optimize, integer, costs, budgets, transfer_smoothness, occupancy_smoothness, admissions_smoothness)
end

function CapacityParams(;optimize=true, costs_setup=false, costs_breakdown=false, costs_convert=false, costs_bedday=false, costs_unitday=false, ordered=false, baseline_always=false, enforce_leadtime=false, capacity_smoothness=0, max_occupancy=1)
    return CapacityParams(optimize, costs_setup, costs_breakdown, costs_convert, costs_bedday, costs_unitday, ordered, baseline_always, enforce_leadtime, capacity_smoothness, max_occupancy)
end

function SolverParams(;timelimit=Inf, verbose=false)
    return SolverParams(timelimit, verbose)
end
