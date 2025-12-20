using CSV
using DataFrames
using Serialization
using Dates
using Distributions
using LinearAlgebra: dot

include("config.jl")

# Default LOS distribution for non-surge patients (shape, scale)
const NONSURGE_LOS = Gamma(0.4872, 8.0417)

alwayszero(x) = 0

"""
    load_capacity(hospitals; patient_type=DEFAULT_PATIENT_TYPE)

Load capacity table for the given patient type and return a vector of
per-hospital capacity levels (beds, priority, etc.).
"""
function load_capacity(hospitals; patient_type=DEFAULT_PATIENT_TYPE)
    capacity_path = joinpath(PROJECT_ROOT, "data", "capacity_$(patient_type).csv")
    capacity_df = DataFrame(CSV.File(capacity_path))

    sort!(capacity_df, [:hospital_id, :priority])
    capacity_data = map(collect(groupby(capacity_df, :hospital_id))) do df
        h = map(eachrow(df)) do row
            location = findfirst(==(row.hospital_id), hospitals)
            capacity_level = findfirst(==(row.capacity_level), CAPACITY_LEVELS)
            (; location, capacity_level, row.beds, row.priority)
        end
        filter(r -> !isnothing(r.capacity_level), h)
    end

    return capacity_data
end

"""
    compute_capacity_params(capacity_data; ...)

Attach cost and timing parameters to each capacity record using the provided
callables. Defaults to zero-valued functions.
"""
function compute_capacity_params(
    capacity_data;
    cost_setup=alwayszero,
    cost_breakdown=alwayszero,
    cost_convert=alwayszero,
    cost_per_unitday=alwayszero,
    cost_per_bedday=alwayszero,
    time_setup=alwayszero,
    time_breakdown=alwayszero,
)
    N = length(capacity_data)
    B = [length(capacity_data[i]) for i in 1:N]

    params = map(1:N) do i
        map(1:B[i]) do b
            r = capacity_data[i][b]
            merge(r, (;
                cost_setup = cost_setup(r),
                cost_breakdown = cost_breakdown(r),
                cost_convert = cost_convert(r),
                cost_per_unitday = cost_per_unitday(r),
                cost_per_bedday = cost_per_bedday(r),
                time_setup = time_setup(r),
                time_breakdown = time_breakdown(r),
            ))
        end
    end

    return params
end

"""
    load_total_capacity(hospitals; patient_type=DEFAULT_PATIENT_TYPE)

Load the total staffed capacity per hospital from the capacity_total CSV file.
Returns `nothing` if the data file is not available.
"""
function load_total_capacity(hospitals; patient_type=DEFAULT_PATIENT_TYPE)
    path = joinpath(PROJECT_ROOT, "data", "capacity_total_$(patient_type).csv")
    if !isfile(path)
        return nothing
    end

    df = DataFrame(CSV.File(path))
    total_capacity = zeros(length(hospitals))
    for (i, h) in enumerate(hospitals)
        idx = findfirst(==(h), df.hospital_id)
        if !isnothing(idx)
            total_capacity[i] = df[idx, :beds]
        end
    end
    return total_capacity
end

"""
    load_nonsurge_data(hospitals, dates_all, surge_arrivals; patient_type=DEFAULT_PATIENT_TYPE)

Load total admissions data and compute non-surge patient census.
Non-surge admissions = total admissions - surge admissions.
Non-surge census is computed using the NONSURGE_LOS distribution.
Returns `nothing` if the data file is not available.
"""
function load_nonsurge_data(hospitals, dates_all, surge_arrivals; patient_type=DEFAULT_PATIENT_TYPE)
    path = joinpath(PROJECT_ROOT, "data", "processed", "daily_admissions.csv")
    if !isfile(path)
        return nothing
    end

    df = DataFrame(CSV.File(path))
    df.date = Date.(df.date)

    N = length(hospitals)
    T = length(dates_all)

    # Build total admissions matrix
    total_admissions = zeros(N, T)
    for row in eachrow(df)
        i = findfirst(==(row.hospital), hospitals)
        t = findfirst(==(row.date), dates_all)
        if !isnothing(i) && !isnothing(t)
            total_admissions[i, t] = row.admissions
        end
    end

    # Non-surge admissions = total - surge
    nonsurge_admissions = max.(0, total_admissions .- surge_arrivals)

    # Compute non-surge census using fixed LOS distribution
    # L[k] = CDF(k-1) = P(LOS <= k-1) = probability of having discharged by day k-1
    # For patient admitted on day s, probability still present on day t is 1 - CDF(t-s)
    # occupancy[t] = sum_{s=1}^{t} admissions[s] * (1 - CDF(t-s))
    #              = sum(admissions[1:t]) - sum_{s=1}^{t} admissions[s] * CDF(t-s)
    L = [cdf(NONSURGE_LOS, t) for t in 0:T]  # L[k+1] = CDF(k) = P(discharged by day k)
    nonsurge_occupancy = zeros(N, T)
    for i in 1:N
        for t in 1:T
            # L[t:-1:1] gives [CDF(t-1), CDF(t-2), ..., CDF(0)]
            # dot with admissions[1:t] gives expected cumulative discharges by day t
            discharges = dot(nonsurge_admissions[i, 1:t], L[t:-1:1])
            nonsurge_occupancy[i, t] = sum(nonsurge_admissions[i, 1:t]) - discharges
        end
    end

    return nonsurge_occupancy
end

"""
    load_data(; start_date, end_date, patient_type=DEFAULT_PATIENT_TYPE)

Load the serialized dataset and trim to the requested date range.
"""
function load_data(; start_date::Date, end_date::Date, patient_type::Symbol=DEFAULT_PATIENT_TYPE)
    rawdata = deserialize(joinpath(PROJECT_ROOT, "data", "data.jlser"))

    hospitals = rawdata.realdata[:meta].hospitals
    dates_all = rawdata.realdata[:meta].date_range

    hospital_names = ["H$i" for i in 1:length(hospitals)]

    dates = filter(d -> start_date <= d <= end_date, dates_all)
    start_date_idx = findfirst(==(start_date), dates_all)
    end_date_idx = findlast(==(end_date), dates_all)
    Topt = collect(start_date_idx:end_date_idx)

    occupancy = rawdata.realdata[patient_type].active
    arrivals = rawdata.realdata[patient_type].admitted

    occupancy = occupancy[:, 1:end_date_idx]
    arrivals = arrivals[:, 1:end_date_idx]

    N = size(arrivals, 1)
    T = size(arrivals, 2)

    capacity = load_capacity(hospitals; patient_type)

    # Load optional non-surge data (returns nothing if data unavailable)
    total_capacity = load_total_capacity(hospitals; patient_type)
    nonsurge_occupancy = load_nonsurge_data(hospitals, dates_all[1:end_date_idx], arrivals; patient_type)

    return (;
        N, T, Topt,
        start_date, end_date,
        hospitals, hospital_names,
        dates_all, dates_opt=dates,
        occupancy, arrivals,
        capacity,
        total_capacity,        # optional: total staffed capacity per hospital
        nonsurge_occupancy,    # optional: non-surge patient census [N, T]
    )
end
