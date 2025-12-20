using CSV
using DataFrames
using Dates
using Distributions
using Statistics

"""
Fit a Gamma distribution to length-of-stay (LOS) data from JHH.

LOS is computed by matching admission and discharge times for each patient.
Uses the same filtering as aggregate_daily_data.jl (excludes Pediatric, Nursery, Emergency Medicine).
"""
function fit_los_distribution(;
    admits_path::String = "data/raw/admits.csv",
    discharges_path::String = "data/raw/discharges.csv",
    hospital::String = "JHH",
)
    println("Loading and processing data for $hospital...")

    # Load data
    admits = load_data(admits_path)
    discharges = load_data(discharges_path)

    println("  Loaded $(nrow(admits)) admits and $(nrow(discharges)) discharges")

    # Filter out excluded services
    admits = filter_services(admits)
    discharges = filter_services(discharges)
    println("  After service filtering: $(nrow(admits)) admits, $(nrow(discharges)) discharges")

    # Extract hospital and filter to target hospital
    admits.hospital = extract_hospital.(admits.department)
    discharges.hospital = extract_hospital.(discharges.department)

    admits = filter(r -> !ismissing(r.hospital) && r.hospital == hospital, admits)
    discharges = filter(r -> !ismissing(r.hospital) && r.hospital == hospital, discharges)
    println("  After hospital filtering: $(nrow(admits)) admits, $(nrow(discharges)) discharges")

    # Compute LOS for matched patients
    los_data = compute_los(admits, discharges)
    println("  Computed LOS for $(length(los_data)) patients")

    # Filter out invalid LOS values
    los_data = filter(x -> x > 0 && isfinite(x), los_data)
    println("  Valid LOS values: $(length(los_data))")

    # Print summary statistics
    println("\nLOS Summary Statistics:")
    println("  Min:    $(round(minimum(los_data), digits=2)) days")
    println("  Max:    $(round(maximum(los_data), digits=2)) days")
    println("  Mean:   $(round(mean(los_data), digits=2)) days")
    println("  Median: $(round(median(los_data), digits=2)) days")
    println("  Std:    $(round(std(los_data), digits=2)) days")

    # Fit Gamma distribution using MLE
    println("\nFitting Gamma distribution...")
    gamma_dist = fit(Gamma, los_data)

    # Extract parameters
    shape = Distributions.shape(gamma_dist)  # α (shape parameter)
    scale = Distributions.scale(gamma_dist)  # θ (scale parameter)

    println("\nGamma Distribution Parameters:")
    println("  Shape (α): $(round(shape, digits=4))")
    println("  Scale (θ): $(round(scale, digits=4))")
    println("  Rate (β = 1/θ): $(round(1/scale, digits=4))")
    println("\n  Distribution: Gamma($(round(shape, digits=4)), $(round(scale, digits=4)))")
    println("  Mean (αθ): $(round(mean(gamma_dist), digits=2)) days")
    println("  Variance (αθ²): $(round(var(gamma_dist), digits=2))")

    # Compute goodness of fit statistics
    println("\nGoodness of Fit:")
    theoretical_quantiles = [0.25, 0.5, 0.75, 0.9, 0.95]
    println("  Quantile comparison (observed vs theoretical):")
    for q in theoretical_quantiles
        obs = quantile(los_data, q)
        theo = quantile(gamma_dist, q)
        println("    $(Int(q*100))th: $(round(obs, digits=2)) vs $(round(theo, digits=2))")
    end

    return gamma_dist, los_data
end

"""
Load admits or discharges CSV, normalizing column names.
"""
function load_data(path::String)
    df = DataFrame(CSV.File(
        path,
        types = Dict(:EncounterEpicCsn => Int),
        silencewarnings = true,
        dateformat = dateformat"yyyy-mm-dd HH:MM:SS",
        missingstring = ["NA", "", "NULL"],
    ))

    # Normalize column names (handle both cases)
    if "startinstant" in names(df)
        rename!(df, :startinstant => :datetime)
    elseif "StartInstant" in names(df)
        rename!(df, :StartInstant => :datetime)
    end

    if "departmentname" in names(df)
        rename!(df, :departmentname => :department)
    elseif "DepartmentName" in names(df)
        rename!(df, :DepartmentName => :department)
    end

    rename!(df,
        :EncounterEpicCsn => :patient_id,
        :HospitalService => :service,
    )

    dropmissing!(df, [:patient_id, :datetime])
    select!(df, [:patient_id, :datetime, :department, :service])

    return df
end

"""
Filter out patients with Pediatric, Nursery, or EMERGENCY MEDICINE in service.
"""
function filter_services(df::DataFrame)
    excluded_patterns = ["Pediatric", "Nursery", "EMERGENCY MEDICINE"]

    filter(df) do row
        ismissing(row.service) && return false
        service = uppercase(row.service)
        for pattern in excluded_patterns
            if occursin(uppercase(pattern), service)
                return false
            end
        end
        return true
    end
end

"""
Extract hospital identifier from department name (first word).
"""
function extract_hospital(department)
    ismissing(department) && return missing
    dept_str = strip(string(department))
    isempty(dept_str) && return missing
    parts = split(dept_str)
    isempty(parts) && return missing
    return String(parts[1])
end

"""
Compute length of stay in days for each patient by matching admits and discharges.
Returns a vector of LOS values in days (as Float64).
"""
function compute_los(admits::DataFrame, discharges::DataFrame)
    # Get first admission datetime per patient
    admits_unique = combine(
        groupby(admits, :patient_id),
        :datetime => minimum => :admit_time,
    )

    # Get last discharge datetime per patient
    discharges_unique = combine(
        groupby(discharges, :patient_id),
        :datetime => maximum => :discharge_time,
    )

    # Inner join to get matched patients
    stays = innerjoin(admits_unique, discharges_unique, on = :patient_id)

    # Compute LOS in days
    los_values = Float64[]
    for row in eachrow(stays)
        admit_time = row.admit_time
        discharge_time = row.discharge_time

        # Skip if discharge is before admission
        if discharge_time < admit_time
            continue
        end

        # Compute LOS in days (as a fraction)
        los_ms = Dates.value(discharge_time - admit_time)  # milliseconds
        los_days = los_ms / (1000 * 60 * 60 * 24)  # convert to days

        push!(los_values, los_days)
    end

    return los_values
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    fit_los_distribution()
end
