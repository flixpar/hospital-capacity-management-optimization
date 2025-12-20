#!/usr/bin/env julia
"""
Compare results from two scenarios within a single experiment folder side by side.

Usage:
    julia scripts/compare_scenarios.jl results-01 base base-transfers [output.html]

This script generates an HTML file showing matching plots from both scenarios
side by side for easy comparison.
"""

using Dates

function find_scenario_pngs(folder::String, scenario::String)
    scenario_dir = joinpath(folder, "png", scenario)
    if !isdir(scenario_dir)
        error("Scenario directory not found: $scenario_dir")
    end

    pngs = Dict{String, String}()
    for f in readdir(scenario_dir)
        if endswith(f, ".png")
            pngs[f] = joinpath(scenario_dir, f)
        end
    end
    return pngs
end

function list_available_scenarios(folder::String)
    png_dir = joinpath(folder, "png")
    if !isdir(png_dir)
        return String[]
    end
    return filter(d -> isdir(joinpath(png_dir, d)), readdir(png_dir))
end

function generate_html(folder::String, scenario1::String, scenario2::String, output_path::String)
    pngs1 = find_scenario_pngs(folder, scenario1)
    pngs2 = find_scenario_pngs(folder, scenario2)

    # Find all unique plot names
    all_plots = sort(unique([keys(pngs1)..., keys(pngs2)...]))

    # Group plots by prefix (e.g., "capacity_", "load_", etc.)
    groups = Dict{String, Vector{String}}()
    for plot in all_plots
        # Extract prefix before first underscore or use "other"
        parts = split(replace(plot, ".png" => ""), "_")
        prefix = length(parts) > 1 ? parts[1] : "other"
        if !haskey(groups, prefix)
            groups[prefix] = String[]
        end
        push!(groups[prefix], plot)
    end

    # Generate HTML
    html = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Scenario Comparison: $(scenario1) vs $(scenario2)</title>
    <style>
        :root {
            --bg-primary: #1a1a2e;
            --bg-secondary: #16213e;
            --bg-card: #0f3460;
            --accent: #e94560;
            --text-primary: #eaeaea;
            --text-secondary: #a0a0a0;
            --border: #2a2a4a;
        }

        * {
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: var(--bg-primary);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem;
        }

        header {
            text-align: center;
            margin-bottom: 3rem;
            padding-bottom: 2rem;
            border-bottom: 1px solid var(--border);
        }

        h1 {
            font-size: 2rem;
            font-weight: 600;
            margin-bottom: 0.5rem;
        }

        .meta {
            color: var(--text-secondary);
            font-size: 0.9rem;
        }

        .folder-info {
            background: var(--bg-card);
            padding: 0.5rem 1rem;
            border-radius: 6px;
            margin-top: 1rem;
            display: inline-block;
            font-size: 0.9rem;
        }

        .folder-labels {
            display: flex;
            justify-content: center;
            gap: 4rem;
            margin-top: 1.5rem;
        }

        .folder-label {
            background: var(--bg-card);
            padding: 0.75rem 2rem;
            border-radius: 8px;
            font-weight: 500;
        }

        .folder-label.left {
            border-left: 4px solid #4ecdc4;
        }

        .folder-label.right {
            border-left: 4px solid #ff6b6b;
        }

        .group {
            margin-bottom: 3rem;
        }

        .group-header {
            background: var(--bg-secondary);
            padding: 1rem 1.5rem;
            border-radius: 8px 8px 0 0;
            border-bottom: 2px solid var(--accent);
            font-size: 1.25rem;
            font-weight: 600;
            text-transform: capitalize;
        }

        .comparison-grid {
            display: grid;
            gap: 1.5rem;
            padding: 1.5rem;
            background: var(--bg-secondary);
            border-radius: 0 0 8px 8px;
        }

        .comparison-row {
            background: var(--bg-card);
            border-radius: 8px;
            overflow: hidden;
        }

        .plot-title {
            padding: 0.75rem 1rem;
            background: rgba(0,0,0,0.2);
            font-weight: 500;
            font-size: 0.95rem;
        }

        .plot-pair {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 1px;
            background: var(--border);
        }

        .plot-container {
            background: white;
            padding: 1rem;
            display: flex;
            align-items: center;
            justify-content: center;
            min-height: 300px;
        }

        .plot-container img {
            max-width: 100%;
            height: auto;
        }

        .plot-container.missing {
            background: var(--bg-secondary);
            color: var(--text-secondary);
            font-style: italic;
        }

        .summary {
            background: var(--bg-secondary);
            padding: 1.5rem;
            border-radius: 8px;
            margin-bottom: 2rem;
        }

        .summary h2 {
            font-size: 1.1rem;
            margin-bottom: 1rem;
            color: var(--accent);
        }

        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
        }

        .stat {
            background: var(--bg-card);
            padding: 1rem;
            border-radius: 6px;
        }

        .stat-value {
            font-size: 1.5rem;
            font-weight: 600;
        }

        .stat-label {
            font-size: 0.85rem;
            color: var(--text-secondary);
        }

        .toc {
            background: var(--bg-secondary);
            padding: 1.5rem;
            border-radius: 8px;
            margin-bottom: 2rem;
        }

        .toc h2 {
            font-size: 1.1rem;
            margin-bottom: 1rem;
        }

        .toc-list {
            display: flex;
            flex-wrap: wrap;
            gap: 0.5rem;
        }

        .toc-list a {
            color: var(--text-primary);
            text-decoration: none;
            background: var(--bg-card);
            padding: 0.5rem 1rem;
            border-radius: 4px;
            transition: background 0.2s;
            text-transform: capitalize;
        }

        .toc-list a:hover {
            background: var(--accent);
        }
    </style>
</head>
<body>
    <header>
        <h1>Scenario Comparison</h1>
        <p class="meta">Generated on $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))</p>
        <div class="folder-info">Results folder: $(basename(abspath(folder)))</div>
        <div class="folder-labels">
            <div class="folder-label left">$(scenario1)</div>
            <div class="folder-label right">$(scenario2)</div>
        </div>
    </header>

    <div class="summary">
        <h2>Summary</h2>
        <div class="summary-grid">
            <div class="stat">
                <div class="stat-value">$(length(pngs1))</div>
                <div class="stat-label">Plots in $(scenario1)</div>
            </div>
            <div class="stat">
                <div class="stat-value">$(length(pngs2))</div>
                <div class="stat-label">Plots in $(scenario2)</div>
            </div>
            <div class="stat">
                <div class="stat-value">$(length(intersect(keys(pngs1), keys(pngs2))))</div>
                <div class="stat-label">Matching plots</div>
            </div>
            <div class="stat">
                <div class="stat-value">$(length(groups))</div>
                <div class="stat-label">Plot groups</div>
            </div>
        </div>
    </div>

    <div class="toc">
        <h2>Plot Groups</h2>
        <div class="toc-list">
"""

    for group_name in sort(collect(keys(groups)))
        html *= """            <a href="#$(group_name)">$(group_name)</a>\n"""
    end

    html *= """
        </div>
    </div>
"""

    # Generate comparison sections for each group
    for group_name in sort(collect(keys(groups)))
        html *= """
    <div class="group" id="$(group_name)">
        <div class="group-header">$(group_name) plots</div>
        <div class="comparison-grid">
"""

        for plot_file in sort(groups[group_name])
            plot_name = replace(plot_file, ".png" => "")

            has1 = haskey(pngs1, plot_file)
            has2 = haskey(pngs2, plot_file)

            # Convert absolute paths to relative paths from output HTML location
            output_dir = dirname(abspath(output_path))
            path1 = has1 ? relpath(pngs1[plot_file], output_dir) : ""
            path2 = has2 ? relpath(pngs2[plot_file], output_dir) : ""

            html *= """
            <div class="comparison-row">
                <div class="plot-title">$(plot_name)</div>
                <div class="plot-pair">
"""

            if has1
                html *= """                    <div class="plot-container"><img src="$(path1)" alt="$(plot_name)"></div>\n"""
            else
                html *= """                    <div class="plot-container missing">Not in $(scenario1)</div>\n"""
            end

            if has2
                html *= """                    <div class="plot-container"><img src="$(path2)" alt="$(plot_name)"></div>\n"""
            else
                html *= """                    <div class="plot-container missing">Not in $(scenario2)</div>\n"""
            end

            html *= """
                </div>
            </div>
"""
        end

        html *= """
        </div>
    </div>
"""
    end

    html *= """
</body>
</html>
"""

    write(output_path, html)
    return output_path
end

function main()
    if length(ARGS) < 3
        println("Usage: julia scripts/compare_scenarios.jl <folder> <scenario1> <scenario2> [output.html]")
        println("")
        println("Arguments:")
        println("  folder      Results folder (e.g., results-01)")
        println("  scenario1   First scenario name (e.g., base)")
        println("  scenario2   Second scenario name (e.g., base-transfers)")
        println("  output.html Output HTML file (default: scenario-comparison.html)")
        println("")

        # If folder is provided, list available scenarios
        if length(ARGS) >= 1 && isdir(ARGS[1])
            scenarios = list_available_scenarios(ARGS[1])
            if !isempty(scenarios)
                println("Available scenarios in $(ARGS[1]):")
                for s in scenarios
                    println("  - $s")
                end
            end
        end
        exit(1)
    end

    folder = ARGS[1]
    scenario1 = ARGS[2]
    scenario2 = ARGS[3]
    output_path = length(ARGS) >= 4 ? ARGS[4] : "scenario-comparison.html"

    # Validate folder exists
    if !isdir(folder)
        error("Folder not found: $folder")
    end

    # Validate scenarios exist
    available = list_available_scenarios(folder)
    if !(scenario1 in available)
        error("Scenario '$scenario1' not found in $folder. Available: $(join(available, ", "))")
    end
    if !(scenario2 in available)
        error("Scenario '$scenario2' not found in $folder. Available: $(join(available, ", "))")
    end

    println("Comparing scenarios in $(abspath(folder)):")
    println("  Left:  $scenario1")
    println("  Right: $scenario2")
    println("")

    output_file = generate_html(folder, scenario1, scenario2, output_path)

    println("Comparison generated: $(abspath(output_file))")
    println("")
    println("Open in browser:")
    println("  open $(output_file)")
end

main()
