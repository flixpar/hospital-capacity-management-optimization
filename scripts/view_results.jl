#!/usr/bin/env julia
"""
View results from a single experiment folder.

Usage:
    julia scripts/view_results.jl results-01 [output.html]

This script generates an HTML file showing all plots from an experiment folder.
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
                rel_path = relpath(joinpath(root, f), png_dir)
                pngs[rel_path] = joinpath(root, f)
            end
        end
    end
    return pngs
end

function generate_html(folder::String, output_path::String)
    pngs = find_pngs(folder)
    all_keys = sort(collect(keys(pngs)))

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

    html = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Experiment Results: $(basename(abspath(folder)))</title>
    <style>
        :root {
            --bg-primary: #1a1a2e;
            --bg-secondary: #16213e;
            --bg-card: #0f3460;
            --accent: #4ecdc4;
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

        .folder-name {
            background: var(--bg-card);
            display: inline-block;
            padding: 0.5rem 1.5rem;
            border-radius: 6px;
            border-left: 4px solid var(--accent);
            margin-top: 1rem;
        }

        .meta {
            color: var(--text-secondary);
            font-size: 0.9rem;
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
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
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
            color: var(--bg-primary);
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

        .plots-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(450px, 1fr));
            gap: 1.5rem;
            padding: 1.5rem;
            background: var(--bg-secondary);
            border-radius: 0 0 8px 8px;
        }

        .plot-card {
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

        .plot-container {
            background: white;
            padding: 1rem;
            display: flex;
            align-items: center;
            justify-content: center;
        }

        .plot-container img {
            max-width: 100%;
            height: auto;
            cursor: pointer;
            transition: transform 0.2s;
        }

        .plot-container img:hover {
            transform: scale(1.02);
        }

        /* Lightbox */
        .lightbox {
            display: none;
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0,0,0,0.9);
            z-index: 1000;
            cursor: pointer;
            align-items: center;
            justify-content: center;
        }

        .lightbox.active {
            display: flex;
        }

        .lightbox img {
            max-width: 95%;
            max-height: 95%;
        }

        .lightbox-close {
            position: absolute;
            top: 1rem;
            right: 1.5rem;
            color: white;
            font-size: 2rem;
            cursor: pointer;
        }
    </style>
</head>
<body>
    <header>
        <h1>Experiment Results</h1>
        <p class="meta">Generated on $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))</p>
        <div class="folder-name">$(basename(abspath(folder)))</div>
    </header>

    <div class="summary">
        <h2>Summary</h2>
        <div class="summary-grid">
            <div class="stat">
                <div class="stat-value">$(length(pngs))</div>
                <div class="stat-label">Total plots</div>
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
        html *= """            <a href="#$(scenario)">$(scenario) ($(length(scenarios[scenario])))</a>\n"""
    end

    html *= """
        </div>
    </div>
"""

    output_dir = dirname(abspath(output_path))

    for scenario in sort(collect(keys(scenarios)))
        html *= """
    <div class="scenario" id="$(scenario)">
        <div class="scenario-header">$(scenario)</div>
        <div class="plots-grid">
"""

        for plot_key in sort(scenarios[scenario])
            plot_name = splitpath(plot_key)[end]
            plot_name = replace(plot_name, ".png" => "")
            rel_path = relpath(pngs[plot_key], output_dir)

            html *= """
            <div class="plot-card">
                <div class="plot-title">$(plot_name)</div>
                <div class="plot-container">
                    <img src="$(rel_path)" alt="$(plot_name)" onclick="openLightbox(this.src)">
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
    <div class="lightbox" id="lightbox" onclick="closeLightbox()">
        <span class="lightbox-close">&times;</span>
        <img id="lightbox-img" src="" alt="">
    </div>

    <script>
        function openLightbox(src) {
            document.getElementById('lightbox-img').src = src;
            document.getElementById('lightbox').classList.add('active');
        }

        function closeLightbox() {
            document.getElementById('lightbox').classList.remove('active');
        }

        document.addEventListener('keydown', function(e) {
            if (e.key === 'Escape') closeLightbox();
        });
    </script>
</body>
</html>
"""

    write(output_path, html)
    return output_path
end

function main()
    if length(ARGS) < 1
        println("Usage: julia scripts/view_results.jl <folder> [output.html]")
        println("")
        println("Arguments:")
        println("  folder      Results folder (e.g., results-01)")
        println("  output.html Output HTML file (default: results.html)")
        exit(1)
    end

    folder = ARGS[1]
    output_path = length(ARGS) >= 2 ? ARGS[2] : "results.html"

    if !isdir(folder)
        error("Folder not found: $folder")
    end

    println("Viewing experiment: $(abspath(folder))")
    println("")

    output_file = generate_html(folder, output_path)

    println("Results page generated: $(abspath(output_file))")
    println("")
    println("Open in browser:")
    println("  open $(output_file)")
end

main()
