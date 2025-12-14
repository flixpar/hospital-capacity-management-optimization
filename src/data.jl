using CSV
using DataFrames
using Serialization
using Dates

include("config.jl")

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

    return (;
        N, T, Topt,
        start_date, end_date,
        hospitals, hospital_names,
        dates_all, dates_opt=dates,
        occupancy, arrivals,
        capacity,
    )
end
