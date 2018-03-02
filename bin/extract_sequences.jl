using NPZ
using PyCall
using MatrixMarket
using ProgressMeter
include("../chaldea/edit_sim.jl")
include("../chaldea/barton_sternberg.jl")
include("../chaldea/LSH.jl")
include("../chaldea/utils/myutils.jl")

function apply_gaussian_filter(mat, σ)
    if σ == 0
        return mat
    end
    normal = Normal(0, σ)
    extent = 3σ
    normal_filter = transpose(pdf(normal, -float(extent):float(extent)))
    rows, cols = findn(mat .> 0)
    filtered_mat = zeros(size(mat)[1], (size(mat)[2]+2extent))
    for (row, col) in zip(rows, cols)
      # はみ出し防止のため，ちょっと大きめに行列を確保している
      # colのズレを修正
      col += extent
      filtered_mat[row, col-extent:col+extent] += normal_filter
    end
    # "みみ"の部分を切り捨てる
    filtered_mat = filtered_mat[:, (extent+1):(end-extent)]
    return(filtered_mat)
end

profiles_dir = ARGS[1]
hosei = parse(Float64, ARGS[2])
sigma = parse(Int64, ARGS[3]) # size of gaussian filter
a_ = 0
d = chomp(readstring(`date '+%Y%m%dT%H%M%S'`))
save_dir = joinpath(profiles_dir, "sequences_sigma"*string(sigma)*"_hosei"*string(hosei)*"_"*"a"*string(a_)*"_"*d)
#save_dir = joinpath(profiles_dir, "sequences_hosei"*string(hosei)*"_"*d)
@show save_dir
if !isdir(save_dir)
    mkdir(save_dir)
end
open(fid->write(fid, "$save_dir"), "/tmp/save_dir", "w")
println("load start")
(profiles, data, correspondence_table_time, binsize, reference_points_list, simmat_data, params, binarray_data) = myutils.ret_data(profiles_dir);
println("load end")
binarray_csc = myutils.ret_coo_matrix_from_data(binarray_data);
num_spikes = zeros(size(binarray_csc)[1])
# todo: 確率の計算をもう一度確認すること
binarray_csc_T = transpose(binarray_csc)
for row in 1:size(binarray_csc)[1]
    num_spikes[row] = sum(binarray_csc_T[:, row])
end
T = size(binarray_csc)[2]
base_line = sum((num_spikes ./ T) .* (num_spikes ./ T))
window = div(params["window(ms)"] , binsize)
@show window
# reference_pointsに対応する部分神経活動の抽出
mats_list = []
seqtonum = div(1000, binsize)
@show reference_points_list
for reference_points in reference_points_list
    mats = []
    for (idx, reference_timing) in enumerate(reference_points)
        rt = ceil(Int, reference_timing * seqtonum)
        rt = rt >= 1 ? rt : 1
        push!(mats, full(binarray_csc[:, (rt):(rt+window)]))
    end
    push!(mats_list, mats)
end

parse_int(x) = parse(Int, x)
idx_list = map(parse_int, keys(profiles))
@show(length(idx_list))
sequences_list = []
reduced_rp_list = []
aligned_mats_list = []
progress = Progress(length(idx_list))
@show mats_list
for (i, idx) in enumerate(idx_list)
    #    profile = apply_gaussian_filter(profiles[key], sigma)
    profile = profiles[string(idx)]
    mats = mats_list[idx]
    reference_points = reference_points_list[idx]
    sequence_rows_list = []
    sequence_cols_list = []
    aligned_mats = []
    for (reference_point, mat) in zip(reference_points, mats)
        #bp, dp, maxrow, maxcol = edit_sim_local_with_gap_with_bp(mat, profile, a=params["a"], base=0)
        bp, dp, maxrow, maxcol = edit_sim_local_with_gap_with_bp(mat, profile, a=a_, base=0)
        mat, _ = find_common(bp, maxrow, maxcol, mat, profile, hosei=hosei)
        (rows, cols) = findn(mat .>= 1)
        cols /= seqtonum
        cols += reference_point
        if length(unique(rows)) >= 5
            push!(aligned_mats, mat)
            push!(sequence_rows_list, rows)
            push!(sequence_cols_list, cols)
        end
    end
    push!(sequences_list, (sequence_rows_list, sequence_cols_list))
    push!(aligned_mats_list, aligned_mats)
    next!(progress)    
end

save_rows_dict = Dict{String, Any}()
save_cols_dict = Dict{String, Any}()
for idx in idx_list
    (sequence_rows_list, sequence_cols_list) = sequences_list[idx]
    idx2 = 1
    for (srl, scl) in zip(sequence_rows_list, sequence_cols_list)
        save_rows_dict[string(idx)*"."*string(idx2)] = srl - 1
        save_cols_dict[string(idx)*"."*string(idx2)] = scl
        idx2 += 1
    end
end

npzwrite(joinpath(save_dir, "sequence_rows.npz"), save_rows_dict)
npzwrite(joinpath(save_dir, "sequence_cols.npz"), save_cols_dict)

