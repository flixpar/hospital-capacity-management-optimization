using DataFrames
using Dates

"""
    compute_metrics(data, results)

Summaries of transfers, occupancy, admissions, and capacity by hospital/day.
"""
function compute_metrics(data, results)
    N, _, T = size(results.transfers)

    metrics_total = (;
        hospitals = length(data.hospitals),
        days = length(data.dates_opt),
        transfers = sum(results.transfers),
        occupancy = sum(results.occupancy),
        admissions = sum(results.admissions),
        allocated_capacity = sum(results.capacity),
    )

    metrics_byhospital = DataFrame([
        (
            hospital = data.hospitals[i],
            hospital_name = data.hospital_names[i],
            transfers_out = sum(results.transfers[i,:,:]),
            transfers_in = sum(results.transfers[:,i,:]),
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
