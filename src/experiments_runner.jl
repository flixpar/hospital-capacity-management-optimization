using DataFrames

include("config.jl")
include("data.jl")
include("los.jl")
include("model_types.jl")
include("model_optimize.jl")
include("metrics.jl")
include("plotting.jl")

"""
    run_experiment(exp_id, capacity_cost_fns, capacity_params, transfer_params, solver_params;
                   data, los_dists, output_path=DEFAULT_OUTPUT_PATH, show=true)

Complete pipeline: load/prepare data, optimize decisions, compute metrics, and
write plots. `data` and `los_dists` can be reused across scenarios.
"""
function run_experiment(
    exp_id,
    capacity_cost_fns,
    capacity_params,
    transfer_params,
    solver_params;
    data,
    los_dists=nothing,
    output_path=DEFAULT_OUTPUT_PATH,
    show=true,
)
    los = isnothing(los_dists) ? estimate_los(data.arrivals, data.occupancy, data.Topt) : los_dists
    capacity_cost_params = compute_capacity_params(data.capacity; capacity_cost_fns...)

    model = optimize_decisions(
        data.arrivals,
        capacity_cost_params,
        los,
        data.Topt,
        capacity_params,
        transfer_params,
        solver_params,
    )

    results = unpack_decisions(model)
    metrics = compute_metrics(data, results)

    plot_capacity(data, results, metrics, show=show, save=true, folder=exp_id, output_path=output_path)
    plot_load(data, results, metrics, 1, show=show, save=true, folder=exp_id, output_path=output_path)
    plot_load(data, results, metrics, :dynamic, show=show, save=true, folder=exp_id, output_path=output_path)

    plot_system_load(data, results, metrics, 1, false, show=false, save=false, folder=exp_id, output_path=output_path)
    plot_system_load(data, results, metrics, :dynamic, false, show=false, save=false, folder=exp_id, output_path=output_path)

    plot_system_load(data, results, metrics, 1, true, show=show, save=true, folder=exp_id, output_path=output_path)
    plot_system_load(data, results, metrics, :dynamic, true, show=show, save=true, folder=exp_id, output_path=output_path)

    plot_surge_timeline(data, results, metrics, show=show, save=true, folder=exp_id, output_path=output_path)

    plot_unit_usage(data, results, metrics, 1, false, show=show, save=true, folder=exp_id, output_path=output_path)
    plot_unit_usage(data, results, metrics, 2, false, show=false, save=true, folder=exp_id, output_path=output_path)
    plot_unit_usage(data, results, metrics, 3, false, show=show, save=true, folder=exp_id, output_path=output_path)
    plot_unit_usage(data, results, metrics, 4, false, show=false, save=true, folder=exp_id, output_path=output_path)
    plot_unit_usage(data, results, metrics, 5, false, show=false, save=true, folder=exp_id, output_path=output_path)

    return (; metrics, results, data)
end

"""
    run_experiment_results(exp_id, capacity_cost_fns, capacity_params, transfer_params, solver_params;
                           data, los_dists, output_path)

Variant that only returns metrics/results without plotting.
"""
function run_experiment_results(
    exp_id,
    capacity_cost_fns,
    capacity_params,
    transfer_params,
    solver_params;
    data,
    los_dists=nothing,
    output_path=DEFAULT_OUTPUT_PATH,
)
    los = isnothing(los_dists) ? estimate_los(data.arrivals, data.occupancy, data.Topt) : los_dists
    capacity_cost_params = compute_capacity_params(data.capacity; capacity_cost_fns...)

    model = optimize_decisions(
        data.arrivals,
        capacity_cost_params,
        los,
        data.Topt,
        capacity_params,
        transfer_params,
        solver_params,
    )

    results = unpack_decisions(model)
    metrics = compute_metrics(data, results)

    return (; metrics, results, data)
end

"""
    default_scenarios(data)

Return a vector of scenario definitions mirroring the notebook experiments.
"""
function default_scenarios(data)
    N, T = data.N, data.T

    return [
        (;
            id = "base",
            capacity_cost_fns = (;
                cost_setup = r -> r.beds / 2,
                cost_per_unitday = r -> ((r.capacity_level / 10) * r.beds) + (0.001 * r.priority),
            ),
            capacity_params = CapacityParams(
                optimize=true,
                costs_setup=true,
                costs_unitday=true,
                ordered=false,
                max_occupancy=0.95,
            ),
            transfer_params = TransferParams(
                N, T;
                optimize=false,
                integer=false,
                costs=0.01,
            ),
            solver_params = SolverParams(
                timelimit=30.0,
                verbose=false,
            ),
            show = false,
        ),
        (;
            id = "base-ordered",
            capacity_cost_fns = (;
                cost_setup = r -> r.beds,
                cost_per_unitday = r -> ((r.capacity_level / 10) * r.beds) + (0.001 * r.priority),
            ),
            capacity_params = CapacityParams(
                optimize=true,
                costs_setup=true,
                costs_unitday=true,
                ordered=true,
                max_occupancy=0.95,
            ),
            transfer_params = TransferParams(
                N, T;
                optimize=false,
                integer=false,
                costs=0.01,
            ),
            solver_params = SolverParams(
                timelimit=30.0,
                verbose=false,
            ),
            show = false,
        ),
        (;
            id = "base-transfers",
            capacity_cost_fns = (;
                cost_setup = r -> r.beds,
                cost_per_unitday = r -> ((r.capacity_level / 10) * r.beds) + (0.001 * r.priority),
            ),
            capacity_params = CapacityParams(
                optimize=true,
                costs_setup=true,
                costs_unitday=true,
                ordered=false,
                max_occupancy=0.95,
            ),
            transfer_params = TransferParams(
                N, T;
                optimize=true,
                integer=true,
                costs=0.01,
            ),
            solver_params = SolverParams(
                timelimit=30.0,
                verbose=false,
            ),
            show = false,
        ),
        (;
            id = "base-impractical",
            capacity_cost_fns = (;
                cost_setup = r -> r.beds,
                cost_per_unitday = r -> ((r.capacity_level / 10) * r.beds) + (0.001 * r.priority),
            ),
            capacity_params = CapacityParams(
                optimize=true,
                costs_setup=false,
                costs_unitday=true,
                ordered=false,
                max_occupancy=0.95,
            ),
            transfer_params = TransferParams(
                N, T;
                optimize=false,
                integer=false,
                costs=0.01,
            ),
            solver_params = SolverParams(
                timelimit=30.0,
                verbose=false,
            ),
            show = false,
        ),
        (;
            id = "base-baselinealways",
            capacity_cost_fns = (;
                cost_setup = r -> (r.capacity_level == 1 ? 0.0 : r.beds),
                cost_per_unitday = r -> ((r.capacity_level / 4) * r.beds) + (0.001 * r.priority),
            ),
            capacity_params = CapacityParams(
                optimize=true,
                costs_setup=true,
                costs_unitday=true,
                baseline_always=true,
                ordered=false,
                max_occupancy=0.95,
            ),
            transfer_params = TransferParams(
                N, T;
                optimize=false,
                integer=false,
                costs=0.1,
            ),
            solver_params = SolverParams(
                timelimit=30.0,
                verbose=false,
            ),
            show = false,
        ),
        (;
            id = "base-baselinealways-transfers",
            capacity_cost_fns = (;
                cost_setup = r -> (r.capacity_level == 1 ? 0.0 : r.beds),
                cost_per_unitday = r -> ((r.capacity_level / 4) * r.beds) + (0.001 * r.priority),
            ),
            capacity_params = CapacityParams(
                optimize=true,
                costs_setup=true,
                costs_unitday=true,
                baseline_always=true,
                ordered=false,
                max_occupancy=0.95,
            ),
            transfer_params = TransferParams(
                N, T;
                optimize=true,
                integer=true,
                costs=0.1,
            ),
            solver_params = SolverParams(
                timelimit=30.0,
                verbose=false,
            ),
            show = false,
        ),
        (;
            id = "base-highsetupcost",
            capacity_cost_fns = (;
                cost_setup = r -> 5 * r.beds,
                cost_per_unitday = r -> ((r.capacity_level / 10) * r.beds) + (0.001 * r.priority),
            ),
            capacity_params = CapacityParams(
                optimize=true,
                costs_setup=true,
                costs_unitday=true,
                ordered=false,
                max_occupancy=0.95,
            ),
            transfer_params = TransferParams(
                N, T;
                optimize=false,
                integer=false,
                costs=0.01,
            ),
            solver_params = SolverParams(
                timelimit=30.0,
                verbose=false,
            ),
            show = false,
        ),
    ]
end

"""
    run_scenarios(data; los_dists, output_path)

Run all default scenarios and return a dictionary keyed by experiment id.
"""
function run_scenarios(data; los_dists=nothing, output_path=DEFAULT_OUTPUT_PATH)
    scenarios = default_scenarios(data)
    results = Dict{String, Any}()

    for s in scenarios
        res = run_experiment(
            s.id,
            s.capacity_cost_fns,
            s.capacity_params,
            s.transfer_params,
            s.solver_params;
            data,
            los_dists,
            output_path,
            show=s.show,
        )
        results[s.id] = res
    end

    return results
end

"""
    surge_level_heuristics(data, results_impractical, tfr_results, base_results)

Replicate heuristic surge-level comparison plots.
"""
function surge_level_heuristics(data, results_impractical, tfr_results, base_results; output_path=DEFAULT_OUTPUT_PATH)
    capacity_by_surge = DataFrame(vcat(data.capacity...))
    capacity_by_surge = combine(
        groupby(capacity_by_surge, [:location, :capacity_level]),
        :beds => sum => :beds,
    )
    C = maximum(capacity_by_surge.capacity_level)
    capacity_by_surge = Dict(
        (i.location, i.capacity_level) => i.beds
        for i in eachrow(capacity_by_surge)
    )
    capacity_by_surge = [get(capacity_by_surge, (h, c), 0) for h in 1:data.N, c in 1:C]
    capacity_by_surge = cumsum(capacity_by_surge, dims=2)

    z = [zeros(Int, 5) capacity_by_surge]
    surge_level_timeline = data.occupancy .> permutedims(stack([z]), (1,3,2))
    surge_level_timeline = diff(surge_level_timeline, dims=3) .== -1
    surge_level_timeline = argmax(surge_level_timeline, dims=3)
    surge_level_timeline = [surge_level_timeline[i,t][3] for i in 1:data.N, t in data.Topt]
    bed_requirements = [capacity_by_surge[i,surge_level_timeline[i,t]] for i in 1:data.N, t in 1:length(data.Topt)]

    surge_timeline_sim = fill(-1, (data.N, length(data.Topt)))
    for i in 1:data.N
        current_state = 1
        for t in 1:length(data.Topt)
            occ = data.occupancy[i, data.Topt[t]]
            min_state = findfirst(>(occ), 0.9 * capacity_by_surge[i,:])
            max_state = findlast(<=(occ), 0.7 * capacity_by_surge[i,:])

            isnothing(min_state) && (min_state = C)
            isnothing(max_state) && (max_state = 1)

            if (occ > 0.9 * capacity_by_surge[i,current_state]) && (current_state < C)
                current_state = min_state
            elseif (occ < 0.7 * capacity_by_surge[i,current_state]) && (current_state > 1)
                if min_state <= max_state
                    current_state = min_state
                end
            end

            surge_timeline_sim[i,t] = current_state
        end
    end
    bed_requirements_sim = [capacity_by_surge[i,surge_timeline_sim[i,t]] for i in 1:data.N, t in 1:length(data.Topt)]

    surge_compare = DataFrame([
        (; scenario = "Optimal bed-level allocation", beddays = sum(data.occupancy[:,data.Topt])),
        (; scenario = "Optimal unit-level allocation", beddays = sum(results_impractical.results.capacity)),
        (; scenario = "Practical unit-level allocation (with transfers)", beddays = sum(tfr_results.results.capacity)),
        (; scenario = "Practical unit-level allocation", beddays = sum(base_results.results.capacity)),
        (; scenario = "Optimal surge-level allocation", beddays = sum(bed_requirements)),
        (; scenario = "Simulated surge-level allocation", beddays = sum(bed_requirements_sim)),
    ])
    surge_compare.pct_diff = surge_compare.beddays / sum(data.occupancy[:,data.Topt])
    surge_compare.pct_diff_alt = surge_compare.beddays / sum(bed_requirements_sim)

    plt = plot(
        reverse(surge_compare),
        x = :beddays,
        y = :scenario,
        Geom.bar(orientation=:horizontal),
        Guide.xlabel("Total Required Dedicated Capacity (Bed-Days)"),
        Guide.ylabel(""),
        style(
            background_color=colorant"white",
            bar_spacing=1mm,
            plot_padding=[0mm, 15mm, 5mm, 0mm];
            fontstyles...,
        ),
    )
    savefig(plt, (32cm, 8cm), output_path, "compare-strategies", "allocation-compare")

    plot_surge_timeline(
        data,
        (;
            capacity_by_surgelevel = capacity_by_surge,
            surge_level = surge_level_timeline,
        ),
        nothing,
        show=true,
        save=true,
        folder="compare-strategies",
        output_path=output_path,
    )

    plot_surge_timeline(
        data,
        (;
            capacity_by_surgelevel = capacity_by_surge,
            surge_level = surge_timeline_sim,
        ),
        nothing,
        show=true,
        save=true,
        folder="compare-strategies",
        output_path=output_path,
    )

    return surge_compare
end
