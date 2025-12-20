function _require(cond::Bool, msg::AbstractString)
    cond || throw(ArgumentError(msg))
    return nothing
end

function _require_finite(name::AbstractString, arr)
    _require(!any(x -> ismissing(x) || !isfinite(x), arr), "$name contains NaN, Inf, or missing values.")
    return nothing
end

function _require_nonnegative(name::AbstractString, arr; allow_inf::Bool=false)
    if allow_inf
        _require(!any(x -> ismissing(x) || isnan(x) || x < 0 || (isinf(x) && x < 0), arr),
                 "$name contains negative, NaN, -Inf, or missing values.")
    else
        _require_finite(name, arr)
        _require(!any(x -> x < 0, arr), "$name contains negative values.")
    end
    return nothing
end

function _require_finite_scalar(name::AbstractString, value; allow_inf::Bool=false)
    _require(!ismissing(value), "$name is missing.")
    if allow_inf
        _require(!isnan(value), "$name is NaN.")
        _require(!(isinf(value) && value < 0), "$name is -Inf.")
    else
        _require(isfinite(value), "$name must be finite.")
    end
    return nothing
end

function validate_optimize_inputs(
    arrivals,
    capacity,
    los,
    Topt,
    capacity_params,
    transfer_params,
    solver_params;
    nonsurge_occupancy=nothing,
    total_capacity=nothing,
)
    _require(ndims(arrivals) == 2, "arrivals must be an N x T matrix.")
    N, T = size(arrivals)
    _require(N > 0 && T > 0, "arrivals must be a non-empty N x T matrix.")
    _require_finite("arrivals", arrivals)
    _require_nonnegative("arrivals", arrivals)

    _require(length(Topt) > 0, "Topt must be non-empty.")
    _require(issorted(Topt), "Topt must be sorted in ascending order.")
    _require(length(unique(Topt)) == length(Topt), "Topt must not contain duplicates.")
    _require(all(1 .<= Topt .<= T), "Topt must be within 1:T.")
    if length(Topt) > 1
        _require(all(diff(Topt) .== 1), "Topt must be a contiguous range (step = 1).")
    end

    _require(length(los) == N, "los must have length N.")
    L = discretize_los(los, N, T)
    _require_finite("los pdf", L)
    _require_nonnegative("los pdf", L)

    _require(length(capacity) == N, "capacity must have length N.")
    for i in 1:N
        _require(!isempty(capacity[i]), "capacity[$i] must contain at least one capacity unit.")
        if capacity_params.baseline_always
            _require(any(r -> haskey(r, :capacity_level) && r.capacity_level == 1, capacity[i]),
                     "capacity[$i] must include a baseline capacity_level == 1 when baseline_always is enabled.")
        end
        for (b, record) in enumerate(capacity[i])
            _require(haskey(record, :beds), "capacity[$i][$b] is missing :beds.")
            _require_finite_scalar("capacity[$i][$b].beds", record.beds)
            _require(record.beds >= 0, "capacity[$i][$b].beds must be nonnegative.")

            if capacity_params.baseline_always
                _require(haskey(record, :capacity_level), "capacity[$i][$b] is missing :capacity_level.")
                _require_finite_scalar("capacity[$i][$b].capacity_level", record.capacity_level)
            end

            if capacity_params.enforce_leadtime
                _require(haskey(record, :time_setup) && haskey(record, :time_breakdown),
                         "capacity[$i][$b] is missing :time_setup or :time_breakdown.")
                _require_finite_scalar("capacity[$i][$b].time_setup", record.time_setup)
                _require_finite_scalar("capacity[$i][$b].time_breakdown", record.time_breakdown)
                _require(record.time_setup >= 0 && record.time_setup == floor(record.time_setup),
                         "capacity[$i][$b].time_setup must be a nonnegative integer.")
                _require(record.time_breakdown >= 0 && record.time_breakdown == floor(record.time_breakdown),
                         "capacity[$i][$b].time_breakdown must be a nonnegative integer.")
            end

            if capacity_params.costs_unitday
                _require(haskey(record, :cost_per_unitday), "capacity[$i][$b] is missing :cost_per_unitday.")
                _require_finite_scalar("capacity[$i][$b].cost_per_unitday", record.cost_per_unitday)
                _require(record.cost_per_unitday >= 0, "capacity[$i][$b].cost_per_unitday must be nonnegative.")
            end
            if capacity_params.costs_setup
                _require(haskey(record, :cost_setup), "capacity[$i][$b] is missing :cost_setup.")
                _require_finite_scalar("capacity[$i][$b].cost_setup", record.cost_setup)
                _require(record.cost_setup >= 0, "capacity[$i][$b].cost_setup must be nonnegative.")
            end
            if capacity_params.costs_breakdown
                _require(haskey(record, :cost_breakdown), "capacity[$i][$b] is missing :cost_breakdown.")
                _require_finite_scalar("capacity[$i][$b].cost_breakdown", record.cost_breakdown)
                _require(record.cost_breakdown >= 0, "capacity[$i][$b].cost_breakdown must be nonnegative.")
            end
            if capacity_params.costs_convert
                _require(haskey(record, :cost_convert), "capacity[$i][$b] is missing :cost_convert.")
                _require_finite_scalar("capacity[$i][$b].cost_convert", record.cost_convert)
                _require(record.cost_convert >= 0, "capacity[$i][$b].cost_convert must be nonnegative.")
            end
        end
    end

    _require_finite_scalar("capacity_params.capacity_smoothness", capacity_params.capacity_smoothness)
    _require(capacity_params.capacity_smoothness >= 0, "capacity_params.capacity_smoothness must be nonnegative.")
    _require_finite_scalar("capacity_params.max_occupancy", capacity_params.max_occupancy)
    _require(capacity_params.max_occupancy > 0, "capacity_params.max_occupancy must be > 0.")
    _require_finite_scalar("capacity_params.shortage_penalty", capacity_params.shortage_penalty)
    _require(capacity_params.shortage_penalty >= 0, "capacity_params.shortage_penalty must be nonnegative.")

    costs = transfer_params.costs
    _require(ndims(costs) == 2 && size(costs, 1) == N && size(costs, 2) == N,
             "transfer_params.costs must be an N x N matrix.")
    _require_finite("transfer_params.costs", costs)

    budgets = transfer_params.budgets
    _require(size(budgets.perhospitalday, 1) == N && size(budgets.perhospitalday, 2) == T,
             "transfer_params.budgets.perhospitalday must be N x T.")
    _require(size(budgets.perhospitalpair, 1) == N && size(budgets.perhospitalpair, 2) == N,
             "transfer_params.budgets.perhospitalpair must be N x N.")
    _require(length(budgets.perday) == T, "transfer_params.budgets.perday must have length T.")
    _require(length(budgets.perhospital) == N, "transfer_params.budgets.perhospital must have length N.")
    _require_nonnegative("transfer_params.budgets.perhospitalday", budgets.perhospitalday; allow_inf=true)
    _require_nonnegative("transfer_params.budgets.perhospitalpair", budgets.perhospitalpair; allow_inf=true)
    _require_nonnegative("transfer_params.budgets.perday", budgets.perday; allow_inf=true)
    _require_nonnegative("transfer_params.budgets.perhospital", budgets.perhospital; allow_inf=true)
    _require_finite_scalar("transfer_params.budgets.total", budgets.total; allow_inf=true)
    _require(budgets.total >= 0 || isinf(budgets.total),
             "transfer_params.budgets.total must be nonnegative.")

    _require_finite_scalar("transfer_params.transfer_smoothness", transfer_params.transfer_smoothness)
    _require(transfer_params.transfer_smoothness >= 0, "transfer_params.transfer_smoothness must be nonnegative.")
    _require_finite_scalar("transfer_params.occupancy_smoothness", transfer_params.occupancy_smoothness)
    _require(transfer_params.occupancy_smoothness >= 0, "transfer_params.occupancy_smoothness must be nonnegative.")
    _require_finite_scalar("transfer_params.admissions_smoothness", transfer_params.admissions_smoothness)
    _require(transfer_params.admissions_smoothness >= 0, "transfer_params.admissions_smoothness must be nonnegative.")

    _require_finite_scalar("solver_params.timelimit", solver_params.timelimit; allow_inf=true)
    _require(solver_params.timelimit >= 0 || isinf(solver_params.timelimit),
             "solver_params.timelimit must be nonnegative.")

    if !isnothing(nonsurge_occupancy)
        _require(ndims(nonsurge_occupancy) == 2, "nonsurge_occupancy must be an N x T matrix.")
        _require(size(nonsurge_occupancy, 1) == N && size(nonsurge_occupancy, 2) == T,
                 "nonsurge_occupancy must be an N x T matrix.")
        _require_nonnegative("nonsurge_occupancy", nonsurge_occupancy)
    end
    if !isnothing(total_capacity)
        _require(length(total_capacity) == N, "total_capacity must have length N.")
        _require_nonnegative("total_capacity", total_capacity)
    end
    if capacity_params.shortage_penalty > 0
        _require(!isnothing(nonsurge_occupancy) && !isnothing(total_capacity),
                 "shortage_penalty > 0 requires nonsurge_occupancy and total_capacity.")
    end

    return nothing
end
