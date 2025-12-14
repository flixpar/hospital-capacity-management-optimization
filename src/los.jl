using Distributions
using LinearAlgebra: norm
using BlackBoxOptim
using ProgressMeter

"""
    discretize_los(los, N, T)

Return LOS probability mass per hospital and day (pdf at 0:T).
"""
function discretize_los(los::Array{<:Distribution,1}, N::Int, T::Int)
    return [pdf(l, t) for l in los, t in 0:T]
end

"""
    estimate_los(admissions, occupancy, Topt)

Estimate length-of-stay distributions hospital by hospital via BlackBoxOptim.
"""
function estimate_los(admissions, occupancy, Topt)
    function estimate_los_single(adm, occ)
        unpack_params(params) = Gamma(params[1], params[2])

        function estimate_occupancy_single(adm, los)
            T = length(adm)
            L = pdf.(los, 0:T)
            dis = @views [dot(adm[1:t], L[t:-1:1]) for t in 1:T]
            occ_est = @views [sum(adm[1:t]) - sum(dis[1:t]) for t in 1:T]
            return occ_est
        end

        function score_func(params)
            los_dist = unpack_params(params)
            occ_est = estimate_occupancy_single(adm, los_dist)[Topt]
            return norm(occ_est - occ, 2)
        end

        param_bounds = [(0.0, 40.0), (0.0, 40.0)]
        timelimit = 10.0

        r = bboptimize(
            score_func,
            SearchRange = param_bounds,
            Method = :adaptive_de_rand_1_bin_radiuslimited,
            TraceMode = :silent,
            MaxTime = timelimit,
            RandomizeRngSeed = false,
            RngSeed = 0,
        )

        best_params = best_candidate(r)
        return unpack_params(best_params)
    end

    los_dists = Distribution[]
    @showprogress for i in 1:size(admissions, 1)
        adm = admissions[i, :]
        occ = occupancy[i, Topt]
        push!(los_dists, estimate_los_single(adm, occ))
    end

    return los_dists
end
