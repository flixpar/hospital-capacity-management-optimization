#!/usr/bin/env julia
"""
Script to visualize non-surge patient data including:
- Total admissions vs surge admissions
- Computed non-surge census
- Total staffed capacity

Outputs an HTML page with embedded plots.

Usage:
    julia scripts/view_nonsurge_data.jl [--start-date 2021-12-01] [--end-date 2022-03-01] [--output report.html]
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CSV
using DataFrames
using Dates
using Distributions
using LinearAlgebra: dot
using Gadfly
using Compose
using Colors, ColorSchemes
import Cairo, Fontconfig
using ArgParse

# Include project config
include(joinpath(@__DIR__, "..", "src", "config.jl"))
include(joinpath(@__DIR__, "..", "src", "data.jl"))

function parse_commandline()
    s = ArgParseSettings(description="View non-surge patient data")
    @add_arg_table! s begin
        "--start-date"
            help = "Start date (YYYY-MM-DD)"
            arg_type = String
            default = "2021-01-01"
        "--end-date"
            help = "End date (YYYY-MM-DD)"
            arg_type = String
            default = "2022-03-01"
        "--output", "-o"
            help = "Output HTML file path"
            arg_type = String
            default = "nonsurge_data_report.html"
    end
    return parse_args(s)
end

"""
    compute_total_admissions(hospitals, dates_all)

Load total admissions from the daily_admissions.csv file.
Returns matrix [N hospitals x T days].
"""
function compute_total_admissions(hospitals, dates_all)
    path = joinpath(PROJECT_ROOT, "data", "processed", "daily_admissions.csv")
    if !isfile(path)
        error("Daily admissions file not found: $path")
    end

    df = DataFrame(CSV.File(path))
    df.date = Date.(df.date)

    N = length(hospitals)
    T = length(dates_all)

    total_admissions = zeros(N, T)
    for row in eachrow(df)
        i = findfirst(==(row.hospital), hospitals)
        t = findfirst(==(row.date), dates_all)
        if !isnothing(i) && !isnothing(t)
            total_admissions[i, t] = row.admissions
        end
    end

    return total_admissions
end

"""
    create_admissions_plot(hospitals, hospital_names, dates, total_admissions, surge_arrivals)

Create a plot comparing total vs surge admissions per hospital.
"""
function create_admissions_plot(hospitals, hospital_names, dates, total_admissions, surge_arrivals)
    date_ticks = date_ticks_for_range(dates[1], dates[end])

    # Build DataFrame for plotting
    df_rows = []
    for i in eachindex(hospitals)
        for (t, d) in enumerate(dates)
            push!(df_rows, (
                hospital = hospital_names[i],
                date = d,
                total = total_admissions[i, t],
                surge = surge_arrivals[i, t],
                nonsurge = max(0, total_admissions[i, t] - surge_arrivals[i, t]),
            ))
        end
    end
    df = DataFrame(df_rows)

    # Create stacked bar plot for each hospital
    plots = []
    for h in hospital_names
        h_df = filter(r -> r.hospital == h, df)

        # Reshape for stacked plot
        stacked_df = vcat(
            DataFrame(date=h_df.date, admissions=h_df.nonsurge, type=fill("Non-Surge", nrow(h_df)), order=fill(1, nrow(h_df))),
            DataFrame(date=h_df.date, admissions=h_df.surge, type=fill("Surge (COVID-19)", nrow(h_df)), order=fill(2, nrow(h_df))),
        )

        plt = plot(
            stacked_df,
            x=:date,
            y=:admissions,
            color=:type,
            Geom.bar(position=:stack),
            Guide.xlabel(""),
            Guide.ylabel("Daily Admissions"),
            Guide.title("$(h) - Admissions Breakdown"),
            Guide.colorkey(title=""),
            Scale.color_discrete_manual(colorant"steelblue", colorant"coral", levels=["Non-Surge", "Surge (COVID-19)"]),
            Guide.xticks(ticks=date_ticks),
            Scale.x_continuous(labels=(d -> Dates.format(d, "u d"))),
            Coord.cartesian(ymin=0),
            style(
                background_color=colorant"white",
                bar_spacing=0mm;
                FONTSTYLES...
            ),
        )
        push!(plots, plt)
    end

    return plots
end

"""
    create_census_plot(hospitals, hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)

Create a plot showing non-surge census, surge census, and total capacity.
"""
function create_census_plot(_hospitals, hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)
    date_ticks = date_ticks_for_range(dates[1], dates[end])

    plots = []
    for (i, h) in enumerate(hospital_names)
        # Build data for this hospital
        df_stacked = vcat(
            DataFrame(
                date=dates,
                census=nonsurge_occupancy[i, :],
                type=fill("Non-Surge", length(dates)),
                order=fill(1, length(dates))
            ),
            DataFrame(
                date=dates,
                census=surge_occupancy[i, :],
                type=fill("Surge (COVID-19)", length(dates)),
                order=fill(2, length(dates))
            ),
        )

        total_census = nonsurge_occupancy[i, :] .+ surge_occupancy[i, :]
        cap = isnothing(total_capacity) ? maximum(total_census) * 1.1 : total_capacity[i]
        y_max = max(cap, maximum(total_census)) * 1.1

        layers = [
            layer(
                df_stacked,
                x=:date,
                y=:census,
                color=:type,
                Geom.bar(position=:stack),
                order=-1,
            ),
        ]

        if !isnothing(total_capacity)
            push!(layers, layer(
                yintercept=[cap],
                Geom.hline(color="black", size=1mm, style=:dash),
                order=100,
            ))
        end

        plt = plot(
            layers...,
            Guide.xlabel(""),
            Guide.ylabel("Hospital Census (Beds)"),
            Guide.title("$(h) - Occupancy Breakdown"),
            Guide.colorkey(title=""),
            Scale.color_discrete_manual(colorant"steelblue", colorant"coral", levels=["Non-Surge", "Surge (COVID-19)"]),
            Coord.cartesian(ymin=0, ymax=y_max),
            Guide.xticks(ticks=date_ticks),
            Scale.x_continuous(labels=(d -> Dates.format(d, "u d"))),
            style(
                background_color=colorant"white",
                bar_spacing=0mm;
                FONTSTYLES...
            ),
        )
        push!(plots, plt)
    end

    return plots
end

"""
    create_system_summary_plot(hospitals, hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)

Create a system-wide summary plot.
"""
function create_system_summary_plot(_, _hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)
    date_ticks = date_ticks_for_range(dates[1], dates[end])

    # Aggregate across hospitals
    total_nonsurge = vec(sum(nonsurge_occupancy, dims=1))
    total_surge = vec(sum(surge_occupancy, dims=1))
    total_cap = isnothing(total_capacity) ? nothing : sum(total_capacity)

    df_stacked = vcat(
        DataFrame(
            date=dates,
            census=total_nonsurge,
            type=fill("Non-Surge", length(dates)),
        ),
        DataFrame(
            date=dates,
            census=total_surge,
            type=fill("Surge (COVID-19)", length(dates)),
        ),
    )

    total_census = total_nonsurge .+ total_surge
    y_max = isnothing(total_cap) ? maximum(total_census) * 1.1 : max(total_cap, maximum(total_census)) * 1.1

    layers = [
        layer(
            df_stacked,
            x=:date,
            y=:census,
            color=:type,
            Geom.bar(position=:stack),
            order=-1,
        ),
    ]

    if !isnothing(total_cap)
        push!(layers, layer(
            yintercept=[total_cap],
            Geom.hline(color="black", size=1mm, style=:dash),
            order=100,
        ))
    end

    plt = plot(
        layers...,
        Guide.xlabel(""),
        Guide.ylabel("System Census (Beds)"),
        Guide.title("System-Wide Occupancy"),
        Guide.colorkey(title=""),
        Scale.color_discrete_manual(colorant"steelblue", colorant"coral", levels=["Non-Surge", "Surge (COVID-19)"]),
        Coord.cartesian(ymin=0, ymax=y_max),
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d"))),
        style(
            background_color=colorant"white",
            bar_spacing=0mm;
            FONTSTYLES...
        ),
    )

    return plt
end

"""
    create_utilization_plot(hospitals, hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)

Create a utilization plot showing capacity usage over time.
"""
function create_utilization_plot(_, hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)
    if isnothing(total_capacity)
        return nothing
    end

    date_ticks = date_ticks_for_range(dates[1], dates[end])

    df_rows = []
    for (i, h) in enumerate(hospital_names)
        for (t, d) in enumerate(dates)
            total = nonsurge_occupancy[i, t] + surge_occupancy[i, t]
            util = total_capacity[i] > 0 ? total / total_capacity[i] : 0.0
            push!(df_rows, (
                hospital = h,
                date = d,
                utilization = util,
            ))
        end
    end
    df = DataFrame(df_rows)

    max_util = maximum(df.utilization)

    plt = plot(
        layer(
            ymin=[0.0], ymax=[1.0],
            Geom.hband, alpha=[0.1],
            style(default_color=colorant"green"),
            order=-2,
        ),
        layer(
            ymin=[1.0], ymax=[max(1.5, max_util * 1.1)],
            Geom.hband, alpha=[0.1],
            style(default_color=colorant"red"),
            order=-2,
        ),
        layer(
            df,
            x=:date,
            y=:utilization,
            color=:hospital,
            Geom.line,
            style(line_width=0.5mm),
        ),
        layer(
            yintercept=[1.0],
            Geom.hline(color="black", size=0.5mm, style=:dash),
        ),
        Guide.xlabel(""),
        Guide.ylabel("Capacity Utilization"),
        Guide.title("Total Capacity Utilization by Hospital"),
        Guide.colorkey(title="Hospital"),
        Coord.cartesian(ymin=0, ymax=max(1.5, max_util * 1.1)),
        Guide.xticks(ticks=date_ticks),
        Scale.x_continuous(labels=(d -> Dates.format(d, "u d"))),
        Scale.y_continuous(labels=x -> "$(round(Int, x*100))%"),
        style(
            background_color=colorant"white";
            FONTSTYLES...
        ),
    )

    return plt
end

"""
    plot_to_svg(plt, width_cm, height_cm)

Convert a Gadfly plot to an SVG string.
"""
function plot_to_svg(plt, width_cm, height_cm)
    io = IOBuffer()
    draw(SVG(io, width_cm * cm, height_cm * cm), plt)
    return String(take!(io))
end

"""
    generate_summary_table(hospitals, hospital_names, dates, total_admissions, surge_arrivals,
                           nonsurge_occupancy, surge_occupancy, total_capacity)

Generate HTML table with summary statistics.
"""
function generate_summary_table(_, hospital_names, _, total_admissions, surge_arrivals,
                                 nonsurge_occupancy, surge_occupancy, total_capacity)
    rows = []

    for (i, h) in enumerate(hospital_names)
        total_adm = sum(total_admissions[i, :])
        surge_adm = sum(surge_arrivals[i, :])
        nonsurge_adm = total_adm - surge_adm

        avg_nonsurge_census = mean(nonsurge_occupancy[i, :])
        avg_surge_census = mean(surge_occupancy[i, :])
        avg_total_census = avg_nonsurge_census + avg_surge_census

        peak_total = maximum(nonsurge_occupancy[i, :] .+ surge_occupancy[i, :])

        cap = isnothing(total_capacity) ? "N/A" : "$(Int(total_capacity[i]))"
        peak_util = if isnothing(total_capacity)
            "N/A"
        else
            "$(round(100 * peak_total / total_capacity[i], digits=1))%"
        end

        push!(rows, """
        <tr>
            <td><strong>$(h)</strong></td>
            <td>$(round(Int, nonsurge_adm))</td>
            <td>$(round(Int, surge_adm))</td>
            <td>$(round(Int, total_adm))</td>
            <td>$(round(avg_nonsurge_census, digits=1))</td>
            <td>$(round(avg_surge_census, digits=1))</td>
            <td>$(round(avg_total_census, digits=1))</td>
            <td>$(round(Int, peak_total))</td>
            <td>$(cap)</td>
            <td>$(peak_util)</td>
        </tr>
        """)
    end

    # System totals
    total_adm_sys = sum(total_admissions)
    surge_adm_sys = sum(surge_arrivals)
    nonsurge_adm_sys = total_adm_sys - surge_adm_sys

    total_nonsurge_census = vec(sum(nonsurge_occupancy, dims=1))
    total_surge_census = vec(sum(surge_occupancy, dims=1))
    total_census = total_nonsurge_census .+ total_surge_census

    cap_sys = isnothing(total_capacity) ? "N/A" : "$(Int(sum(total_capacity)))"
    peak_util_sys = if isnothing(total_capacity)
        "N/A"
    else
        "$(round(100 * maximum(total_census) / sum(total_capacity), digits=1))%"
    end

    push!(rows, """
    <tr class="system-total">
        <td><strong>SYSTEM</strong></td>
        <td>$(round(Int, nonsurge_adm_sys))</td>
        <td>$(round(Int, surge_adm_sys))</td>
        <td>$(round(Int, total_adm_sys))</td>
        <td>$(round(mean(total_nonsurge_census), digits=1))</td>
        <td>$(round(mean(total_surge_census), digits=1))</td>
        <td>$(round(mean(total_census), digits=1))</td>
        <td>$(round(Int, maximum(total_census)))</td>
        <td>$(cap_sys)</td>
        <td>$(peak_util_sys)</td>
    </tr>
    """)

    return """
    <table class="summary-table">
        <thead>
            <tr>
                <th>Hospital</th>
                <th>Non-Surge Adm.</th>
                <th>Surge Adm.</th>
                <th>Total Adm.</th>
                <th>Avg Non-Surge Census</th>
                <th>Avg Surge Census</th>
                <th>Avg Total Census</th>
                <th>Peak Census</th>
                <th>Staffed Capacity</th>
                <th>Peak Utilization</th>
            </tr>
        </thead>
        <tbody>
            $(join(rows, "\n"))
        </tbody>
    </table>
    """
end

"""
    generate_html_report(; output_path, plots, summary_table, start_date, end_date, los_info)

Generate the complete HTML report with all plots and tables.
"""
function generate_html_report(; output_path, admissions_plots, census_plots, system_plot,
                               utilization_plot, summary_table, start_date, end_date, los_info)
    # Convert plots to SVG
    admissions_svgs = [plot_to_svg(p, 22, 8) for p in admissions_plots]
    census_svgs = [plot_to_svg(p, 22, 10) for p in census_plots]
    system_svg = plot_to_svg(system_plot, 24, 10)
    util_svg = isnothing(utilization_plot) ? "" : plot_to_svg(utilization_plot, 24, 10)

    html = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Non-Surge Patient Data Report</title>
        <style>
            body {
                font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
                max-width: 1400px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f5f5f5;
            }
            h1, h2, h3 {
                color: #333;
            }
            .header {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 30px;
                border-radius: 10px;
                margin-bottom: 30px;
            }
            .header h1 {
                color: white;
                margin: 0 0 10px 0;
            }
            .header p {
                margin: 5px 0;
                opacity: 0.9;
            }
            .section {
                background: white;
                border-radius: 10px;
                padding: 20px;
                margin-bottom: 20px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            .section h2 {
                border-bottom: 2px solid #667eea;
                padding-bottom: 10px;
                margin-top: 0;
            }
            .plot-container {
                overflow-x: auto;
                margin: 20px 0;
            }
            .plot-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(600px, 1fr));
                gap: 20px;
            }
            .plot-item {
                background: #fafafa;
                padding: 15px;
                border-radius: 8px;
                border: 1px solid #eee;
            }
            .summary-table {
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                font-size: 14px;
            }
            .summary-table th, .summary-table td {
                padding: 12px 15px;
                text-align: right;
                border-bottom: 1px solid #ddd;
            }
            .summary-table th {
                background-color: #667eea;
                color: white;
                font-weight: 600;
            }
            .summary-table td:first-child, .summary-table th:first-child {
                text-align: left;
            }
            .summary-table tr:hover {
                background-color: #f5f5f5;
            }
            .summary-table .system-total {
                background-color: #e8f4f8;
                font-weight: 600;
            }
            .info-box {
                background: #e8f4f8;
                border-left: 4px solid #667eea;
                padding: 15px;
                margin: 15px 0;
                border-radius: 0 8px 8px 0;
            }
            .info-box h4 {
                margin: 0 0 10px 0;
                color: #333;
            }
            .info-box p {
                margin: 5px 0;
            }
            code {
                background: #f0f0f0;
                padding: 2px 6px;
                border-radius: 3px;
                font-family: 'Menlo', 'Monaco', monospace;
                font-size: 0.9em;
            }
            .timestamp {
                color: #666;
                font-size: 0.9em;
            }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Non-Surge Patient Data Report</h1>
            <p><strong>Date Range:</strong> $(Dates.format(start_date, "U d, Y")) to $(Dates.format(end_date, "U d, Y"))</p>
            <p class="timestamp">Generated: $(Dates.format(now(), "U d, Y HH:MM"))</p>
        </div>

        <div class="section">
            <h2>Summary Statistics</h2>
            $(summary_table)

            <div class="info-box">
                <h4>Length of Stay Distribution (Non-Surge Patients)</h4>
                <p>$(los_info)</p>
                <p>This distribution is used to compute non-surge patient census from admissions data.</p>
            </div>
        </div>

        <div class="section">
            <h2>System-Wide Overview</h2>
            <div class="plot-container">
                $(system_svg)
            </div>
            $(util_svg == "" ? "" : "<h3>Capacity Utilization</h3><div class=\"plot-container\">$(util_svg)</div>")
        </div>

        <div class="section">
            <h2>Hospital Census Breakdown</h2>
            <p>Census computed using admission data convolved with length-of-stay distribution.
               Dashed line shows total staffed capacity.</p>
            <div class="plot-grid">
                $(join(["<div class=\"plot-item\">$(svg)</div>" for svg in census_svgs], "\n"))
            </div>
        </div>

        <div class="section">
            <h2>Daily Admissions Breakdown</h2>
            <p>Total daily admissions decomposed into surge (COVID-19) and non-surge patients.</p>
            <div class="plot-grid">
                $(join(["<div class=\"plot-item\">$(svg)</div>" for svg in admissions_svgs], "\n"))
            </div>
        </div>

        <div class="section">
            <h2>Data Sources</h2>
            <div class="info-box">
                <h4>Input Files</h4>
                <p><code>data/processed/daily_admissions.csv</code> - Total hospital admissions</p>
                <p><code>data/data.jlser</code> - Surge patient arrivals and occupancy</p>
                <p><code>data/capacity_total_icu.csv</code> - Total staffed capacity per hospital</p>
            </div>
            <div class="info-box">
                <h4>Computation Method</h4>
                <p><strong>Non-surge admissions</strong> = Total admissions - Surge admissions</p>
                <p><strong>Non-surge census</strong> = Convolution of non-surge admissions with LOS distribution</p>
                <p>Census at time t: <code>sum(admissions[1:t] - discharges[1:t])</code></p>
                <p>Discharges use CDF differences to model minimum 1-day LOS.</p>
            </div>
        </div>
    </body>
    </html>
    """

    open(output_path, "w") do f
        write(f, html)
    end

    println("Report generated: $(abspath(output_path))")
end

function main()
    args = parse_commandline()

    start_date = Date(args["start-date"])
    end_date = Date(args["end-date"])
    output_path = args["output"]

    println("Loading data from $(start_date) to $(end_date)...")

    # Load the main data
    data = load_data(; start_date, end_date)

    hospitals = data.hospitals
    hospital_names = data.hospital_names
    N = data.N

    # Get date range that matches the loaded data
    dates = data.dates_all[1:data.T]

    # Load total admissions
    total_admissions = compute_total_admissions(hospitals, dates)

    # Surge arrivals from loaded data
    surge_arrivals = data.arrivals

    # Non-surge occupancy (computed by load_data using the Gamma LOS distribution)
    if isnothing(data.nonsurge_occupancy)
        error("Non-surge occupancy data not available. Check that data/processed/daily_admissions.csv exists.")
    end
    nonsurge_occupancy = data.nonsurge_occupancy

    # Surge occupancy from the data
    surge_occupancy = data.occupancy

    # Total capacity
    total_capacity = data.total_capacity

    println("Data loaded:")
    println("  Hospitals: $(N)")
    println("  Days: $(length(dates))")
    println("  Date range: $(dates[1]) to $(dates[end])")
    println("  Total capacity available: $(isnothing(total_capacity) ? "No" : "Yes")")

    # LOS distribution info
    los_shape = 0.4872
    los_scale = 8.0417
    los_mean = los_shape * los_scale
    los_info = "Gamma(shape=$(los_shape), scale=$(los_scale)), mean = $(round(los_mean, digits=1)) days"

    println("\nGenerating plots...")

    # Create plots
    admissions_plots = create_admissions_plot(hospitals, hospital_names, dates, total_admissions, surge_arrivals)
    census_plots = create_census_plot(hospitals, hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)
    system_plot = create_system_summary_plot(hospitals, hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)
    utilization_plot = create_utilization_plot(hospitals, hospital_names, dates, nonsurge_occupancy, surge_occupancy, total_capacity)

    # Generate summary table
    summary_table = generate_summary_table(hospitals, hospital_names, dates, total_admissions, surge_arrivals,
                                           nonsurge_occupancy, surge_occupancy, total_capacity)

    println("Generating HTML report...")

    # Generate HTML report
    generate_html_report(;
        output_path,
        admissions_plots,
        census_plots,
        system_plot,
        utilization_plot,
        summary_table,
        start_date = dates[1],
        end_date = dates[end],
        los_info,
    )

    println("\nDone!")
end

main()
