using DataFrames
using Dates

"""
    compute_metrics(data, results)

Summaries of transfers, occupancy, admissions, and capacity by hospital/day.
Includes shortage metrics when non-surge data is available.
"""
function compute_metrics(data, results)
    N, _, T = size(results.transfers)

    # Check if shortage data is available
    has_shortage = !isnothing(results.shortage)
    has_nonsurge = haskey(data, :nonsurge_occupancy) && !isnothing(data.nonsurge_occupancy)
    has_total_capacity = haskey(data, :total_capacity) && !isnothing(data.total_capacity)

    # Base metrics
    metrics_total = (;
        hospitals = length(data.hospitals),
        days = length(data.dates_opt),
        transfers = sum(results.transfers),
        occupancy = sum(results.occupancy),
        admissions = sum(results.admissions),
        allocated_capacity = sum(results.capacity),
        # Shortage metrics (only when data available)
        total_shortage = has_shortage ? sum(results.shortage) : nothing,
        peak_shortage = has_shortage ? maximum(results.shortage) : nothing,
        shortage_days = has_shortage ? sum(results.shortage .> 0.01) : nothing,
        nonsurge_occupancy = has_nonsurge ? sum(data.nonsurge_occupancy[:, data.Topt]) : nothing,
    )

    metrics_byhospital = DataFrame([
        (
            hospital = data.hospitals[i],
            hospital_name = data.hospital_names[i],
            transfers_out = sum(results.transfers[i,:,:]),
            transfers_in = sum(results.transfers[:,i,:]),
            # Shortage metrics by hospital
            total_shortage = has_shortage ? sum(results.shortage[i,:]) : nothing,
            peak_shortage = has_shortage ? maximum(results.shortage[i,:]) : nothing,
            shortage_days = has_shortage ? sum(results.shortage[i,:] .> 0.01) : nothing,
            total_capacity = has_total_capacity ? data.total_capacity[i] : nothing,
        )
        for i in 1:N
    ])

    metrics_bydate = DataFrame([
        (
            date = data.start_date + Day(t-1),
            occupancy = sum(results.occupancy[:,t]),
            admissions = sum(results.admissions[:,t]),
            transfers = sum(results.transfers[:,:,t]),
        )
        for t in 1:T
    ])

    metrics_byhospitalday = DataFrame([
        (
            hospital = data.hospitals[i],
            hospital_name = data.hospital_names[i],
            date = data.dates_opt[t],
            admissions = results.admissions[i,t],
            transfers_out = sum(results.transfers[i,:,t]),
            transfers_in = sum(results.transfers[:,i,t]),
            occupancy = results.occupancy[i,t],
            capacity = results.capacity[i,t],
        )
        for i in 1:N, t in 1:T
    ][:])

    return (;
        total = metrics_total,
        byhospital = metrics_byhospital,
        bydate = metrics_bydate,
        byhospitalday = metrics_byhospitalday,
    )
end
