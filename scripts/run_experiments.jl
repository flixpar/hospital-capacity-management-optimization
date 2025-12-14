using CSV
using DataFrames
using Distributions
using Serialization
using Dates
using LinearAlgebra
using JuMP
using Gurobi
using MathOptInterface
using BlackBoxOptim
using ProgressMeter
using Gadfly
using Compose
using Colors, ColorSchemes

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "experiments_runner.jl"))

# Core study settings (matches the original notebook defaults)
start_date = Date(2021, 12, 15)
end_date = Date(2022, 2, 15)
patient_type = DEFAULT_PATIENT_TYPE
output_path = DEFAULT_OUTPUT_PATH

data = load_data(; start_date, end_date, patient_type)
los_dists = estimate_los(data.arrivals, data.occupancy, data.Topt)

println("Running default scenarios...")
scenario_results = run_scenarios(data; los_dists, output_path)

base_results = scenario_results["base"]
tfr_results = scenario_results["base-transfers"]
impractical_results = scenario_results["base-impractical"]

println("Total number of transfers (base-transfers): $(sum(tfr_results.results.transfers))")

println("Building surge-level comparison plots...")
surge_compare = surge_level_heuristics(
    data,
    impractical_results,
    tfr_results,
    base_results;
    output_path=output_path,
)

println("Done. Outputs saved under $(output_path).")
