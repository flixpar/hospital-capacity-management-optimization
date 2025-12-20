#!/usr/bin/env julia
"""
Compare results from two experiment folders side by side.

Usage:
    julia scripts/compare_results.jl results-01 results-02 [output.html]

This script generates an HTML file showing matching plots from both experiment
folders side by side for easy comparison.
"""

using Dates

function find_pngs(folder::String)
    png_dir = joinpath(folder, "png")
    if !isdir(png_dir)
        error("No 'png' directory found in $folder")
    end

    pngs = Dict{String, String}()
    for (root, dirs, files) in walkdir(png_dir)
        for f in files
            if endswith(f, ".png")
                # Get relative path from png folder: scenario/filename.png
                rel_path = relpath(joinpath(root, f), png_dir)
                pngs[rel_path] = joinpath(root, f)
            end
        end
    end
    return pngs
end

function generate_html(folder1::String, folder2::String, output_path::String)
    pngs1 = find_pngs(folder1)
    pngs2 = find_pngs(folder2)

    # Find all unique keys and organize by scenario
    all_keys = sort(unique([keys(pngs1)..., keys(pngs2)...]))

    # Group by scenario
    scenarios = Dict{String, Vector{String}}()
    for k in all_keys
        parts = splitpath(k)
        scenario = length(parts) > 1 ? parts[1] : "root"
        if !haskey(scenarios, scenario)
            scenarios[scenario] = String[]
        end
        push!(scenarios[scenario], k)
    end

    # Generate HTML
    html = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Experiment Comparison: $(basename(folder1)) vs $(basename(folder2))</title>
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

        .scenario {
            margin-bottom: 3rem;
        }

        .scenario-header {
            background: var(--bg-secondary);
            padding: 1rem 1.5rem;
            border-radius: 8px 8px 0 0;
            border-bottom: 2px solid var(--accent);
            font-size: 1.25rem;
            font-weight: 600;
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
        }

        .toc-list a:hover {
            background: var(--accent);
        }
    </style>
</head>
<body>
    <header>
        <h1>Experiment Comparison</h1>
        <p class="meta">Generated on $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))</p>
        <div class="folder-labels">
            <div class="folder-label left">$(basename(abspath(folder1)))</div>
            <div class="folder-label right">$(basename(abspath(folder2)))</div>
        </div>
    </header>

    <div class="summary">
        <h2>Summary</h2>
        <div class="summary-grid">
            <div class="stat">
                <div class="stat-value">$(length(pngs1))</div>
                <div class="stat-label">Plots in $(basename(folder1))</div>
            </div>
            <div class="stat">
                <div class="stat-value">$(length(pngs2))</div>
                <div class="stat-label">Plots in $(basename(folder2))</div>
            </div>
            <div class="stat">
                <div class="stat-value">$(length(intersect(keys(pngs1), keys(pngs2))))</div>
                <div class="stat-label">Matching plots</div>
            </div>
            <div class="stat">
                <div class="stat-value">$(length(scenarios))</div>
                <div class="stat-label">Scenarios</div>
            </div>
        </div>
    </div>

    <div class="toc">
        <h2>Scenarios</h2>
        <div class="toc-list">
"""

    for scenario in sort(collect(keys(scenarios)))
        html *= """            <a href="#$(scenario)">$(scenario)</a>\n"""
    end

    html *= """
        </div>
    </div>
"""

    # Generate comparison sections for each scenario
    for scenario in sort(collect(keys(scenarios)))
        html *= """
    <div class="scenario" id="$(scenario)">
        <div class="scenario-header">$(scenario)</div>
        <div class="comparison-grid">
"""

        for plot_key in sort(scenarios[scenario])
            plot_name = splitpath(plot_key)[end]
            plot_name = replace(plot_name, ".png" => "")

            has1 = haskey(pngs1, plot_key)
            has2 = haskey(pngs2, plot_key)

            # Convert absolute paths to relative paths from output HTML location
            output_dir = dirname(abspath(output_path))
            path1 = has1 ? relpath(pngs1[plot_key], output_dir) : ""
            path2 = has2 ? relpath(pngs2[plot_key], output_dir) : ""

            html *= """
            <div class="comparison-row">
                <div class="plot-title">$(plot_name)</div>
                <div class="plot-pair">
"""

            if has1
                html *= """                    <div class="plot-container"><img src="$(path1)" alt="$(plot_name)"></div>\n"""
            else
                html *= """                    <div class="plot-container missing">Not in $(basename(folder1))</div>\n"""
            end

            if has2
                html *= """                    <div class="plot-container"><img src="$(path2)" alt="$(plot_name)"></div>\n"""
            else
                html *= """                    <div class="plot-container missing">Not in $(basename(folder2))</div>\n"""
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
    if length(ARGS) < 2
        println("Usage: julia scripts/compare_results.jl <folder1> <folder2> [output.html]")
        println("")
        println("Arguments:")
        println("  folder1     First results folder (e.g., results-01)")
        println("  folder2     Second results folder (e.g., results-02)")
        println("  output.html Output HTML file (default: comparison.html)")
        exit(1)
    end

    folder1 = ARGS[1]
    folder2 = ARGS[2]
    output_path = length(ARGS) >= 3 ? ARGS[3] : "comparison.html"

    # Validate folders exist
    if !isdir(folder1)
        error("Folder not found: $folder1")
    end
    if !isdir(folder2)
        error("Folder not found: $folder2")
    end

    println("Comparing experiments:")
    println("  Left:  $(abspath(folder1))")
    println("  Right: $(abspath(folder2))")
    println("")

    output_file = generate_html(folder1, folder2, output_path)

    println("Comparison generated: $(abspath(output_file))")
    println("")
    println("Open in browser:")
    println("  open $(output_file)")
end

main()
