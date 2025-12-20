using Gadfly
using Compose
using Colors, ColorSchemes
using Dates
using DataFrames
import Cairo, Fontconfig

include("config.jl")

defaultfont = DEFAULT_FONT
const fontstyles = FONTSTYLES

function savefig(plt, size, output_path, folder, filename)
    @assert !isnothing(folder)

    output_path_eps = joinpath(output_path, "eps", folder, filename * ".eps")
    output_path_pdf = joinpath(output_path, "pdf", folder, filename * ".pdf")
    output_path_png = joinpath(output_path, "png", folder, filename * ".png")

    mkpath(dirname(output_path_eps))
    mkpath(dirname(output_path_pdf))
    mkpath(dirname(output_path_png))

    plt |> PS(output_path_eps , size...)
    plt |> PDF(output_path_pdf, size...)
    plt |> PNG(output_path_png, size..., dpi=250)

    return
end

function plot_capacity(data, results, metrics; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    date_ticks = date_ticks_for_range(data.start_date, data.end_date)
    plt = plot(
        metrics.byhospitalday,
        x = :date,
        y = :capacity,
        color = :hospital_name,
        Geom.line,
        Guide.xlabel(""),
        Guide.ylabel("COVID-19 ICU Capacity"),
        Guide.colorkey(title="Hospital"),
        Coord.cartesian(ymin=0, xmin=date_ticks[1], xmax=date_ticks[end]),
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        Scale.y_continuous(format=:plain),
        layer(
            yintercept=[0],
            xintercept=[date_ticks[1]],
            Geom.hline(color="lightgrey", size=0.6mm),
            Geom.vline(color="lightgrey", size=0.6mm),
        ),
        style(
            line_width=0.55mm,
            background_color=colorant"white",
            plot_padding=[5mm, 10mm, 5mm, 5mm];
            fontstyles...,
        ),
    )

    plt_size = (20cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "capacity")
    end

    return plt, plt_size
end

function plot_load(data, results, metrics, capacity_level; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    date_ticks = date_ticks_for_range(data.start_date, data.end_date)
    metrics_df = select(metrics.byhospitalday, :hospital, :hospital_name, :date, :occupancy, :capacity)

    if capacity_level == :dynamic
        metrics_df.load_l = metrics_df.occupancy ./ metrics_df.capacity
    else
        capacity_l_data = [filter(r -> r.capacity_level <= capacity_level, h) for h in data.capacity]
        capacity_l = [sum([r.beds for r in h]) for h in capacity_l_data]

        hosp_idx_lookup = Dict(zip(data.hospitals, 1:data.N))
        metrics_df.hospital_idx = map(h -> hosp_idx_lookup[h], metrics_df.hospital)
        metrics_df.capacity_l = map(h -> capacity_l[h], metrics_df.hospital_idx)
        metrics_df.load_l = metrics_df.occupancy ./ metrics_df.capacity_l
    end

    ylabel = (capacity_level == 1) ? "Baseline Capacity Utilization" : "Occupancy"
    load_plot_max = max(2.0, maximum(metrics_df.load_l))

    plt = plot(
        layer(
            metrics_df,
            x = :date,
            y = :load_l,
            color = :hospital_name,
            Geom.line,
        ),
        Guide.xlabel(""),
        Guide.ylabel(ylabel),
        layer(ymin=[0.0], ymax=[1.0], Geom.hband, alpha=[0.1], style(default_color=colorant"green"), order=-1),
        layer(ymin=[1.0], ymax=[load_plot_max], Geom.hband, alpha=[0.1], style(default_color=colorant"red"), order=-1),
        Coord.cartesian(ymin=0, ymax=load_plot_max, xmin=date_ticks[1], xmax=date_ticks[end]),
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        Scale.y_continuous(labels=x -> "$(round(Int, x*100))%"),
        style(
            background_color=colorant"white",
            plot_padding=[5mm, 10mm, 5mm, 5mm],
            key_position=:none;
            fontstyles...,
        ),
    )

    plt_size = (18cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "load_$(capacity_level)")
    end

    return plt, plt_size
end

function plot_system_load(data, results, metrics, capacity_level, show_all; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    date_ticks = date_ticks_for_range(data.start_date, data.end_date)

    system_df = combine(
        groupby(metrics.byhospitalday, :date),
        :occupancy => sum => :occupancy,
        :capacity => sum => :capacity,
    )

    if capacity_level == :dynamic
        system_df.load_l = system_df.occupancy ./ system_df.capacity
    else
        capacity_l_data = [filter(r -> r.capacity_level <= capacity_level, h) for h in data.capacity]
        capacity_l = [sum([r.beds for r in h]) for h in capacity_l_data]
        capacity_l = sum(capacity_l)
        system_df.capacity_l = fill(capacity_l, nrow(system_df))
        system_df.load_l = system_df.occupancy ./ capacity_l
    end

    metrics_df = select(metrics.byhospitalday, :hospital, :hospital_name, :date, :occupancy, :capacity)

    if capacity_level == :dynamic
        metrics_df.load_l = metrics_df.occupancy ./ metrics_df.capacity
    else
        capacity_l_data = [filter(r -> r.capacity_level <= capacity_level, h) for h in data.capacity]
        capacity_l = [sum([r.beds for r in h]) for h in capacity_l_data]

        hosp_idx_lookup = Dict(zip(data.hospitals, 1:data.N))
        metrics_df.hospital_idx = map(h -> hosp_idx_lookup[h], metrics_df.hospital)
        metrics_df.capacity_l = map(h -> capacity_l[h], metrics_df.hospital_idx)
        metrics_df.load_l = metrics_df.occupancy ./ metrics_df.capacity_l
    end

    sort!(metrics_df, [:hospital_name, :date])
    hospitals = sort(unique(metrics_df.hospital_name))

    ylabel = (capacity_level == 1) ? "Baseline Capacity Utilization" : "Occupancy"
    load_plot_max = max(1.25, maximum(metrics_df.load_l))

    indiv_layer = if show_all
        layer(
            metrics_df,
            x = :date,
            y = :load_l,
            color = :hospital_name,
            Geom.line,
            style(line_width=0.25mm),
            order=-100,
        )
    else
        load_plot_max = 1.2
        Guide.yticks(ticks=0:0.2:1.2)
    end

    colorkey = if show_all
        labels = ["System"; hospitals]
        keycolors = [colorant"navy"; Scale.default_discrete_colors(length(hospitals))]
        Guide.manual_color_key("", labels, keycolors)
    else
        style()
    end

    plt = plot(
        layer(
            system_df,
            x = :date,
            y = :load_l,
            Geom.point,
            Geom.line,
            style(
                line_width=0.7mm,
                default_color=colorant"navy",
            ),
            order=100,
        ),
        layer(
            yintercept = [1.0],
            Geom.hline(color="red", size=0.8mm),
            order=0,
        ),
        layer(
            yintercept = [0.95],
            Geom.hline(color="gold", size=0.8mm),
            order=0,
        ),
        indiv_layer,
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        Scale.y_continuous(format=:plain),
        Scale.y_continuous(labels=x -> "$(round(Int, x*100))%"),
        Guide.xlabel(""),
        Guide.ylabel(ylabel),
        colorkey,
        Coord.cartesian(ymin=0, ymax=load_plot_max, xmin=date_ticks[1], xmax=date_ticks[end]),
        style(
            background_color=colorant"white",
            colorkey_swatch_shape=:square,
            plot_padding=[5mm, 10mm, 5mm, 5mm],
            key_position=:none;
            fontstyles...,
        ),
    )

    plt_size = (20cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "system_load_$(capacity_level)")
    end

    return plt, plt_size
end

function plot_unit_usage(data, results, metrics, h_idx, binary; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    date_ticks = date_ticks_for_range(data.start_date, data.end_date)
    B = length(data.capacity[h_idx])

    if h_idx == 3
        B = length([h for h in data.capacity[h_idx] if h.capacity_level < 4])
    end

    usage = results.capacity_unit_allocated[h_idx,:]
    usage_df = allcombinations(DataFrame, :unit => 1:B, :day => 1:length(usage))
    usage_df.date = map(d -> data.dates_opt[d], usage_df.day)

    if binary
        usage_df.allocated = map(r -> max(usage[r.day][r.unit], 0), eachrow(usage_df))
        sort!(usage_df, :allocated)
    else
        usage_df.allocated = map(r -> (usage[r.day][r.unit] > 0) ? data.capacity[h_idx][r.unit].beds : 0, eachrow(usage_df))
    end

    if binary
        colorscale = Scale.color_discrete_manual(colorant"gray95", colorant"deepskyblue")
        colorkey = Guide.colorkey(title="", labels=["Not Used", "Used"])
    else
        max_z = maximum(b.beds for b in data.capacity[h_idx])
        cmap = ColorScheme(range(colorant"white", colorant"deepskyblue", length=100))
        cmap_fn = z -> (z == 0) ? colorant"gray95" : cmap[0.75*z + 0.25]
        colorscale = Scale.color_continuous(colormap=cmap_fn)
        colorkey = Guide.colorkey(title="Beds")
    end

    yticks = if B < 10
        1:B
    elseif B < 50
        z = collect(0:5:B)
        z[1] = 1
        z
    else
        z = collect(0:10:B)
        z[1] = 1
        z
    end

    levels = unique(b.capacity_level for b in data.capacity[h_idx])
    level_inds = [findlast(b -> b.capacity_level == c, data.capacity[h_idx]) for c in levels]
    level_inds = [i for i in level_inds if !isnothing(i)]
    level_inds = isempty(level_inds) ? [-1] : level_inds
    level_names = [CAPACITY_NAMES[l] for l in levels]

    plt = plot(
        usage_df,
        x = :date,
        y = :unit,
        color = :allocated,
        Geom.rectbin,
        layer(
            yintercept = (level_inds .+ 0.5),
            Geom.hline(color="gray", size=3px),
            x = [data.dates_opt[end] for _ in level_inds],
            y = [y+0.45-(0.03B) for y in level_inds],
            label = level_names,
            Geom.label(position=:left),
            order=100,
        ),
        Guide.xlabel(""),
        Guide.ylabel("Unit ID (H$(h_idx))"),
        Guide.title("H$(h_idx)"),
        colorkey,
        colorscale,
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        Coord.cartesian(ymin=0.5, ymax=B+0.5),
        Guide.yticks(ticks=yticks),
        style(
            background_color=colorant"white",
            plot_padding=[5mm, 5mm, 5mm, 5mm];
            fontstyles...,
            key_label_font_size=14px,
        ),
    )

    plt_size = (24cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "unit_usage_$(h_idx)")
    end

    return plt, plt_size
end

function plot_surge_timeline(data, results, metrics; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    date_ticks = date_ticks_for_range(data.start_date, data.end_date)

    if !haskey(results, :surge_level)
        max_alloc_timeline = map(z -> findlast(>(0),z), results.capacity_unit_allocated)
        capacity_levels = [[b.capacity_level for b in h] for h in data.capacity]
        capacity_timeline = [capacity_levels[i][max_alloc_timeline[i,t]] for i in 1:data.N, t in 1:length(data.Topt)]
    else
        capacity_timeline = results.surge_level
    end

    timeline_df = DataFrame(
        hospital = repeat(data.hospital_names, 1, length(data.Topt))[:],
        date = permutedims(repeat(data.dates_opt, 1, data.N), (2,1))[:],
        date_offset = permutedims(repeat(data.dates_opt, 1, data.N), (2,1))[:] .+ Day(1),
        capacity_level = capacity_timeline[:],
        H0 = fill(0, data.N * length(data.Topt)),
        H1 = fill(1, data.N * length(data.Topt)),
    )

    C = maximum(timeline_df.capacity_level)
    capacity_colors = [colorant"rgb(39,123,70)", colorant"yellow", colorant"gold", colorant"red", colorant"purple", colorant"black"][1:C]

    plt = plot(
        timeline_df,
        xmin = :date,
        xmax = :date_offset,
        ymin = :H0,
        ymax = :H1,
        ygroup = :hospital,
        color = :capacity_level,
        Geom.subplot_grid(
            layer(Geom.rect),
            Guide.yticks(ticks=[0]),
            Guide.xticks(ticks=date_ticks),
        ),
        Guide.xlabel(""),
        Guide.ylabel("Hospital"),
        Guide.colorkey(title="Surge Level", labels=CAPACITY_NAMES[1:C]),
        Scale.color_discrete_manual(capacity_colors..., levels=1:C),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        Scale.y_continuous(labels=(x -> "")),
        style(
            background_color=colorant"white",
            plot_padding=[8mm, 5mm, 5mm, 5mm];
            fontstyles...,
            key_label_font_size=14px,
        ),
    )

    plt_size = (24cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "surge_timeline")
    end

    return plt, plt_size
end

function plot_surge_timeline_alt(data, results)
    date_ticks = date_ticks_for_range(data.start_date, data.end_date)
    max_alloc_timeline = map(z -> findlast(>(0),z), results.capacity_unit_allocated)
    capacity_levels = [[b.capacity_level for b in h] for h in data.capacity]
    capacity_timeline = [capacity_levels[i][max_alloc_timeline[i,t]] for i in 1:data.N, t in 1:length(data.Topt)]

    timeline_df = DataFrame(
        hospital = repeat(data.hospital_names, 1, length(data.Topt))[:],
        date = permutedims(repeat(data.dates_opt, 1, data.N), (2,1))[:],
        capacity_level = capacity_timeline[:],
    )

    C = maximum(timeline_df.capacity_level)
    capacity_colors = [colorant"rgb(39,123,70)", colorant"yellow", colorant"gold", colorant"red", colorant"purple", colorant"black"][1:C]

    plt = plot(
        timeline_df,
        x = :date,
        y = :hospital,
        color = :capacity_level,
        Geom.rectbin,
        Guide.xlabel(""),
        Guide.ylabel("Hospital"),
        Guide.colorkey(title="Surge Level", labels=CAPACITY_NAMES[1:C]),
        Scale.color_discrete_manual(capacity_colors..., levels=1:C),
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        style(
            background_color=colorant"white",
            plot_padding=[5mm, 5mm, 5mm, 5mm];
            fontstyles...,
            key_label_font_size=14px,
        ),
    )

    # plt |> SVG(24cm, 10cm)
    return plt
end

# =============================================================================
# Non-Surge Patient Shortage Plots
# =============================================================================

"""
    plot_shortage_timeline(data, results, metrics; ...)

Show shortage evolution over time for each hospital.
Returns (nothing, nothing) if shortage data is not available.
"""
function plot_shortage_timeline(data, results, metrics; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    if isnothing(results.shortage)
        return nothing, nothing
    end

    date_ticks = date_ticks_for_range(data.start_date, data.end_date)

    shortage_df = DataFrame([
        (
            hospital_name = data.hospital_names[i],
            date = data.dates_opt[t],
            shortage = results.shortage[i,t],
        )
        for i in 1:data.N, t in 1:length(data.Topt)
    ][:])

    max_shortage = max(1.0, maximum(shortage_df.shortage))

    plt = plot(
        shortage_df,
        x = :date,
        y = :shortage,
        color = :hospital_name,
        Geom.line,
        Guide.xlabel(""),
        Guide.ylabel("Non-Surge Patient Shortage (Beds)"),
        Guide.colorkey(title="Hospital"),
        Coord.cartesian(ymin=0, ymax=max_shortage * 1.1, xmin=date_ticks[1], xmax=date_ticks[end]),
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        Scale.y_continuous(format=:plain),
        layer(
            yintercept=[0],
            xintercept=[date_ticks[1]],
            Geom.hline(color="lightgrey", size=0.6mm),
            Geom.vline(color="lightgrey", size=0.6mm),
        ),
        style(
            line_width=0.55mm,
            background_color=colorant"white",
            plot_padding=[5mm, 10mm, 5mm, 5mm];
            fontstyles...,
        ),
    )

    plt_size = (20cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "shortage_timeline")
    end

    return plt, plt_size
end

"""
    plot_occupancy_breakdown(data, results, metrics, h_idx; ...)

Show surge + non-surge patient occupancy as stacked bars with total capacity threshold line.
Returns (nothing, nothing) if non-surge or total capacity data is not available.
"""
function plot_occupancy_breakdown(data, results, metrics, h_idx; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    has_nonsurge = haskey(data, :nonsurge_occupancy) && !isnothing(data.nonsurge_occupancy)
    has_total_capacity = haskey(data, :total_capacity) && !isnothing(data.total_capacity)

    if !has_nonsurge || !has_total_capacity
        return nothing, nothing
    end

    date_ticks = date_ticks_for_range(data.start_date, data.end_date)

    # Build stacked data: non-surge first (bottom), then surge (top)
    T_opt = length(data.Topt)

    nonsurge_df = DataFrame([
        (
            date = data.dates_opt[t],
            patient_type = "Non-Surge",
            census = data.nonsurge_occupancy[h_idx, data.Topt[t]],
            order = 1,
        )
        for t in 1:T_opt
    ])

    surge_df = DataFrame([
        (
            date = data.dates_opt[t],
            patient_type = "Surge (COVID-19)",
            census = results.occupancy[h_idx, t],
            order = 2,
        )
        for t in 1:T_opt
    ])

    breakdown_df = vcat(nonsurge_df, surge_df)
    sort!(breakdown_df, [:date, :order])

    total_cap = data.total_capacity[h_idx]
    max_census = maximum(nonsurge_df.census .+ surge_df.census)
    y_max = max(total_cap, max_census) * 1.1

    plt = plot(
        layer(
            breakdown_df,
            x = :date,
            y = :census,
            color = :patient_type,
            Geom.bar(position=:stack),
            order = -1,
        ),
        layer(
            yintercept = [total_cap],
            Geom.hline(color="black", size=1mm, style=:dash),
            order = 100,
        ),
        Guide.xlabel(""),
        Guide.ylabel("Hospital Census (Beds)"),
        Guide.title("$(data.hospital_names[h_idx]) - Occupancy Breakdown"),
        Guide.colorkey(title="Patient Type"),
        Scale.color_discrete_manual(colorant"steelblue", colorant"coral", levels=["Non-Surge", "Surge (COVID-19)"]),
        Coord.cartesian(ymin=0, ymax=y_max, xmin=date_ticks[1], xmax=date_ticks[end]),
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        Scale.y_continuous(format=:plain),
        Guide.manual_color_key("", ["Total Staffed Capacity"], [colorant"black"]),
        style(
            background_color=colorant"white",
            bar_spacing=0mm,
            plot_padding=[5mm, 10mm, 5mm, 5mm];
            fontstyles...,
        ),
    )

    plt_size = (20cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "occupancy_breakdown_$(h_idx)")
    end

    return plt, plt_size
end

"""
    plot_system_shortage(data, results, metrics; ...)

Show aggregate shortage across the entire hospital system.
Returns (nothing, nothing) if shortage data is not available.
"""
function plot_system_shortage(data, results, metrics; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    if isnothing(results.shortage)
        return nothing, nothing
    end

    date_ticks = date_ticks_for_range(data.start_date, data.end_date)

    system_shortage = vec(sum(results.shortage, dims=1))

    shortage_df = DataFrame(
        date = data.dates_opt,
        shortage = system_shortage,
        zero = fill(0.0, length(system_shortage)),
    )

    max_shortage = max(1.0, maximum(system_shortage))

    plt = plot(
        layer(
            shortage_df,
            x = :date,
            ymin = :zero,
            ymax = :shortage,
            Geom.ribbon,
            style(default_color=colorant"rgba(180,50,50,0.3)"),
            order = -1,
        ),
        layer(
            shortage_df,
            x = :date,
            y = :shortage,
            Geom.line,
            Geom.point,
            style(
                line_width=0.7mm,
                default_color=colorant"darkred",
                point_size=2pt,
            ),
            order = 100,
        ),
        layer(
            yintercept=[0],
            Geom.hline(color="lightgrey", size=0.6mm),
        ),
        Guide.xlabel(""),
        Guide.ylabel("System-Wide Shortage (Beds)"),
        Coord.cartesian(ymin=0, ymax=max_shortage * 1.1, xmin=date_ticks[1], xmax=date_ticks[end]),
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d, Y"))),
        Scale.y_continuous(format=:plain),
        style(
            background_color=colorant"white",
            plot_padding=[5mm, 10mm, 5mm, 5mm];
            fontstyles...,
        ),
    )

    plt_size = (20cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "system_shortage")
    end

    return plt, plt_size
end

"""
    plot_scenario_comparison(results_base, results_penalty; ...)

Compare key metrics between two scenarios (e.g., with vs without shortage penalty).
"""
function plot_scenario_comparison(results_base, results_penalty; show=false, save=false, folder=nothing, output_path=DEFAULT_OUTPUT_PATH)
    # Extract metrics from both scenarios
    base_shortage = isnothing(results_base.metrics.total.total_shortage) ? 0.0 : results_base.metrics.total.total_shortage
    base_peak = isnothing(results_base.metrics.total.peak_shortage) ? 0.0 : results_base.metrics.total.peak_shortage
    base_days = isnothing(results_base.metrics.total.shortage_days) ? 0 : results_base.metrics.total.shortage_days
    base_capacity = sum(results_base.results.capacity)

    penalty_shortage = isnothing(results_penalty.metrics.total.total_shortage) ? 0.0 : results_penalty.metrics.total.total_shortage
    penalty_peak = isnothing(results_penalty.metrics.total.peak_shortage) ? 0.0 : results_penalty.metrics.total.peak_shortage
    penalty_days = isnothing(results_penalty.metrics.total.shortage_days) ? 0 : results_penalty.metrics.total.shortage_days
    penalty_capacity = sum(results_penalty.results.capacity)

    comparison_df = DataFrame([
        (scenario = "Base (No Penalty)", metric = "Total Shortage", value = base_shortage),
        (scenario = "Base (No Penalty)", metric = "Peak Shortage", value = base_peak),
        (scenario = "Base (No Penalty)", metric = "Shortage Days", value = Float64(base_days)),
        (scenario = "Base (No Penalty)", metric = "Capacity (Bed-Days)", value = base_capacity / 1000),  # Scale down for visibility
        (scenario = "With Shortage Penalty", metric = "Total Shortage", value = penalty_shortage),
        (scenario = "With Shortage Penalty", metric = "Peak Shortage", value = penalty_peak),
        (scenario = "With Shortage Penalty", metric = "Shortage Days", value = Float64(penalty_days)),
        (scenario = "With Shortage Penalty", metric = "Capacity (Bed-Days)", value = penalty_capacity / 1000),
    ])

    plt = plot(
        comparison_df,
        x = :metric,
        y = :value,
        color = :scenario,
        Geom.bar(position=:dodge),
        Guide.xlabel(""),
        Guide.ylabel("Value (Capacity in thousands)"),
        Guide.title("Scenario Comparison: Impact of Shortage Penalty"),
        Guide.colorkey(title="Scenario"),
        Scale.color_discrete_manual(colorant"steelblue", colorant"coral", levels=["Base (No Penalty)", "With Shortage Penalty"]),
        Coord.cartesian(ymin=0),
        Scale.y_continuous(format=:plain),
        style(
            background_color=colorant"white",
            bar_spacing=2mm,
            plot_padding=[5mm, 10mm, 5mm, 5mm];
            fontstyles...,
        ),
    )

    plt_size = (22cm, 10cm)

    if show
        plt |> SVG(plt_size...)
    end
    if save
        savefig(plt, plt_size, output_path, folder, "scenario_comparison")
    end

    return plt, plt_size
end
