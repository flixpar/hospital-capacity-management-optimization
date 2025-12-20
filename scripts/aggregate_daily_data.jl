using CSV
using DataFrames
using Dates

# Main hospitals to include in output
const MAIN_HOSPITALS = ["JHH", "BMC", "HCGH", "SMH", "SH"]

# ICU-related hospital services (adult only, pediatric/neonatal already filtered out)
const ICU_SERVICES = [
    "Intensive Care",
    "Critical Care",
    "Medical ICU",
    "Surgical ICU",
    "Cardiac ICU",
    "Neuro ICU",
    "Burn",
    "Trauma",
]

# JHH ICU departments (JHH codes ICU by location, not service)
const JHH_ICU_DEPARTMENTS = [
    "JHH ZAYED 10E",      # Medical ICU (MICU)
    "JHH ZAYED 9E",       # Surgical ICU (SICU)
    "JHH ZAYED 5E",       # Cardiac ICU (CICU)
    "JHH ZAYED 12W BRU",  # Brain Rescue Unit / Neuro ICU
    "JHH WBG OVFLW PACU ICU",  # Overflow ICU
]

"""
Transform patient-level admits and discharges data into daily hospital-level aggregates.

Outputs:
- data/processed/daily_admissions.csv: Daily admission counts per hospital
- data/processed/daily_census.csv: Daily census (patients present) per hospital
"""
function aggregate_daily_data(;
    admits_path::String = "data/raw/admits.csv",
    discharges_path::String = "data/raw/discharges.csv",
    output_dir::String = "data/processed",
    hospitals::Vector{String} = MAIN_HOSPITALS,
    icu_only::Bool = false,
)
    mkpath(output_dir)

    println("Loading admits data...")
    admits = load_admits(admits_path)
    println("  Loaded $(nrow(admits)) admit records")

    println("Loading discharges data...")
    discharges = load_discharges(discharges_path)
    println("  Loaded $(nrow(discharges)) discharge records")

    println("Filtering out Pediatric, Nursery, and Emergency Medicine patients...")
    admits = filter_services(admits)
    discharges = filter_services(discharges)
    println("  Admits after filtering: $(nrow(admits))")
    println("  Discharges after filtering: $(nrow(discharges))")

    println("Extracting hospital identifiers...")
    admits.hospital = extract_hospital.(admits.department)
    discharges.hospital = extract_hospital.(discharges.department)

    # Remove rows where hospital extraction failed
    admits = filter(r -> !ismissing(r.hospital), admits)
    discharges = filter(r -> !ismissing(r.hospital), discharges)
    println("  Admits with valid hospital: $(nrow(admits))")
    println("  Discharges with valid hospital: $(nrow(discharges))")

    # Filter to main hospitals only
    println("Filtering to main hospitals: $(join(hospitals, ", "))...")
    admits = filter(r -> r.hospital in hospitals, admits)
    discharges = filter(r -> r.hospital in hospitals, discharges)
    println("  Admits after hospital filter: $(nrow(admits))")
    println("  Discharges after hospital filter: $(nrow(discharges))")

    # Filter to ICU services only if requested
    if icu_only
        println("Filtering to ICU services only...")
        admits = filter_icu_services(admits)
        discharges = filter_icu_services(discharges)
        println("  Admits after ICU filter: $(nrow(admits))")
        println("  Discharges after ICU filter: $(nrow(discharges))")
    end

    println("Computing daily admissions...")
    daily_admissions = compute_daily_admissions(admits)
    println("  Generated $(nrow(daily_admissions)) daily admission records")

    println("Computing daily census...")
    daily_census = compute_daily_census(admits, discharges)
    println("  Generated $(nrow(daily_census)) daily census records")

    # Save outputs
    admissions_path = joinpath(output_dir, "daily_admissions.csv")
    census_path = joinpath(output_dir, "daily_census.csv")

    daily_admissions |> CSV.write(admissions_path)
    daily_census |> CSV.write(census_path)

    println("\nOutputs saved to:")
    println("  $admissions_path")
    println("  $census_path")

    return daily_admissions, daily_census
end

"""
Load admits data, normalizing column names.
"""
function load_admits(path::String)
    df = DataFrame(CSV.File(
        path,
        types = Dict(:EncounterEpicCsn => Int),
        silencewarnings = true,
        dateformat = dateformat"yyyy-mm-dd HH:MM:SS",
        missingstring = ["NA", "", "NULL"],
    ))

    # Normalize column names (handle both cases seen in data)
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

    # Drop rows with missing patient_id or datetime
    dropmissing!(df, [:patient_id, :datetime])

    # Extract date from datetime
    df.date = Date.(df.datetime)

    select!(df, [:patient_id, :datetime, :date, :department, :service])
    return df
end

"""
Load discharges data, normalizing column names.
"""
function load_discharges(path::String)
    df = DataFrame(CSV.File(
        path,
        types = Dict(:EncounterEpicCsn => Int),
        silencewarnings = true,
        dateformat = dateformat"yyyy-mm-dd HH:MM:SS",
        missingstring = ["NA", "", "NULL"],
    ))

    # Normalize column names
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

    # Drop rows with missing patient_id or datetime
    dropmissing!(df, [:patient_id, :datetime])

    # Extract date from datetime
    df.date = Date.(df.datetime)

    select!(df, [:patient_id, :datetime, :date, :department, :service])
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
Filter to ICU-related patients only.
For JHH: uses department-based filtering (ICU coded by location)
For other hospitals: uses service-based filtering
"""
function filter_icu_services(df::DataFrame)
    filter(df) do row
        hospital = get(row, :hospital, missing)

        # JHH: filter by department name (ICU is coded by location, not service)
        if hospital == "JHH"
            ismissing(row.department) && return false
            return row.department in JHH_ICU_DEPARTMENTS
        end

        # Other hospitals: filter by service name
        ismissing(row.service) && return false
        return row.service in ICU_SERVICES
    end
end

"""
Extract hospital identifier from department name (first word).
Returns missing if department is missing or empty.
"""
function extract_hospital(department)
    ismissing(department) && return missing
    dept_str = strip(string(department))
    isempty(dept_str) && return missing

    # First word is the hospital identifier
    parts = split(dept_str)
    isempty(parts) && return missing
    return String(parts[1])
end

"""
Compute daily admission counts per hospital.
Returns a long-format DataFrame with columns: date, hospital, admissions
"""
function compute_daily_admissions(admits::DataFrame)
    # Group by date and hospital, count admissions
    daily = combine(
        groupby(admits, [:date, :hospital]),
        nrow => :admissions
    )

    sort!(daily, [:hospital, :date])
    return daily
end

"""
Compute daily census per hospital by matching admits with discharges.
Census on a given day = patients admitted on or before that day who have not yet discharged.
Only includes patients with matching discharge records for accurate census.
Returns a long-format DataFrame with columns: date, hospital, census
"""
function compute_daily_census(admits::DataFrame, discharges::DataFrame)
    # Build patient stays by joining admits and discharges
    # Use admission hospital as the hospital of record

    # Get first admission per patient (in case of multiple records)
    admits_unique = combine(
        groupby(admits, :patient_id),
        :date => minimum => :admit_date,
        :hospital => first => :hospital,
    )

    # Get discharge date per patient
    discharges_unique = combine(
        groupby(discharges, :patient_id),
        :date => maximum => :discharge_date,
    )

    # Inner join to only include patients with both admit and discharge
    # This ensures we have complete stay information for accurate census
    stays = innerjoin(admits_unique, discharges_unique, on = :patient_id)

    # Filter out stays with invalid dates
    filter!(r -> !ismissing(r.admit_date) && !ismissing(r.discharge_date), stays)

    # Get date range
    max_date = maximum(stays.discharge_date)

    # Remove invalid stays where discharge is before admission
    filter!(r -> r.discharge_date >= r.admit_date, stays)

    println("  Built $(nrow(stays)) valid patient stays")

    # Get date range and hospitals
    min_date = minimum(stays.admit_date)
    hospitals = sort(unique(stays.hospital))

    println("  Date range: $min_date to $max_date")
    println("  Hospitals: $(join(hospitals, ", "))")

    # Compute census for each day and hospital
    # This is O(days * hospitals * patients) but more memory efficient
    all_dates = min_date:Day(1):max_date
    census_records = Vector{NamedTuple{(:date, :hospital, :census), Tuple{Date, String, Int}}}()

    for hospital in hospitals
        hospital_stays = filter(r -> r.hospital == hospital, stays)

        for d in all_dates
            # Count patients present on this day
            # Present if: admit_date <= d <= discharge_date
            census = count(r -> r.admit_date <= d <= r.discharge_date, eachrow(hospital_stays))
            push!(census_records, (date = d, hospital = hospital, census = census))
        end
    end

    daily_census = DataFrame(census_records)
    sort!(daily_census, [:hospital, :date])

    return daily_census
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    aggregate_daily_data(icu_only = true)
end
