using NPZ
using MatrixMarket
using ProgressMeter

@everywhere include("../chaldea/edit_sim.jl")
@everywhere include("../chaldea/barton_sternberg.jl")
@everywhere include("../chaldea/LSH.jl")
# redirect stdout for logging
# originalSTDOUT = STDOUT
# (outRead, outWrite) = redirect_stdout()
# println("parallel environment", nworkers()," ", nprocs())
println("calc start")
binarray_file = ARGS[1]
source_dir = dirname(binarray_file)
a = parse(Float64, ARGS[2])   # unit less
window = parse(Int64, ARGS[3]) # (ms)
slide = ARGS[4] # "all" or "half"
min_len_sequence = parse(Int64, ARGS[5]) # unit less
exhaustive = parse(Int64, ARGS[6]) # 0: non-exhaustive, 1: exhaustive
println("loading binarray...")
binarray_data = npzread(binarray_file)
println("binarray loaded")
row = binarray_data["row"] + 1
col = binarray_data["col"] + 1
data = binarray_data["data"]
binsize = binarray_data["binsize(ms)"]
# window(ms)がbinarray上で何column分に対応するのかを格納する．edit_simの計算に必要
window_num = div(window, binsize)
if slide == "half"
    slidewidth = Int64(div(window_num, 2))
else
    slidewidth = window_num
end
duration = binarray_data["duration(ms)"]
row = [Int64(r) for r in row] + 1;
col = [Int64(c) for c in col] + 1;
binarray_csc = sparse(row, col, data)
@show binsize, size(binarray_csc)

println("binarray_file: ", binarray_file)
println("Parameters", "window:", window, "min_len_sequence: ", min_len_sequence, "slidewidth: ", slidewidth)
d = chomp(readstring(`date '+%Y%m%dT%H%M%S'`))
save_dir = joinpath(source_dir, "simmat_window_"*string(window)*"a_"*string(a)*"min_len_"*string(min_len_sequence)*"_"*d)
if !isdir(save_dir)
    mkdir(save_dir)
end
open(fid->write(fid, "$save_dir"), "/tmp/save_dir", "w")
println("savedir: ", save_dir)

binarray_csc_T = transpose(binarray_csc)
# シーケンスが存在する場合のJaccard係数の見積もりと，それを基にしたLSHのBand設計
(nrows, ncols) = size(binarray_csc)
# シーケンスが存在する場合のJaccard係数の見積もりと，それを基にしたLSHのBand設計
num_spikes = zeros(nrows)
# todo: 確率の計算をもう一度確認すること
for row in 1:size(binarray_csc)[1]
    num_spikes[row] = sum(binarray_csc_T[:, row])
end
num_spikes
prob_fire = zeros(size(binarray_csc)[1])
T = ncols
for row in 1:size(binarray_csc)[1]
    # 1ビンに発火が存在する確率
    # (num_spikes[row] / T)  // * binsize?
    # window 個数のビンの少なくとも１つに発火の存在する確率
    # 1-(window個数のビン全てに発火が存在しない確率)
    prob_fire[row] = 1 - (1 - (num_spikes[row] / duration) * binsize)^window_num
end
# # 1つのビンに着目したとき，発火しているニューロンの個数の期待値
NI = sum(prob_fire)
# # 2つのビンに着目したとき，発火している共通ニューロンの個数の期待値
NII = sum(prob_fire .^ 2)
@show NI, NII
# シーケンスが存在しないと仮定した場合のJaccard係数の見積もり
Jaccard1 = (NII) / (2*NI-NII)
# min_len_sequence個のニューロンからなるシーケンスが存在したと仮定した場合の見積もり
Jaccard2 = (NII+min_len_sequence) / (2*NI-NII+min_len_sequence)
@show Jaccard1, Jaccard2
function determine_band(s1, s2, threshold1=.1, threshold2=.7)
    # 1: 雑音
    # 2: シーケンス
    # r: band width
    # b: num band
    println("try...")
    @show threshold1, threshold2
    for r in 1:50
        for b in 1:50
            prob1 = cal_prob(s1, b, r)
            prob2 = cal_prob(s2, b, r)
            if prob1 - threshold1 < 0 && prob2 - threshold2 > 0
                return(r, b)
            end
        end
    end
    println("cannot determine bandwidth and numband!!! Try different variable...")
    # close stdout
    # close(outWrite)
    # data = bytestring(readavailable(outRead))
    # close(outRead)
    # redirect_stdout(originalSTDOUT)
    #
    # open(joinpath(save_dir, "log.txt")) do wf
    #   write(wf, data)
    # end
    threshold1 += 0.05
    determine_band(s1, s2, threshold1, threshold2)
end
#upperbound, lowerbound = 0.5, 0.7
times = 1:slidewidth:(size(binarray_csc)[2] - 2window_num)
npzwrite(joinpath(save_dir, "times.npz"), (collect(times)-1)*binsize)
@eval @everywhere window_num = $window_num
@eval @everywhere slidewidth = $slidewidth
@eval @everywhere min_len_sequence = $min_len_sequence
@eval @everywhere binarray_csc = $binarray_csc
if exhaustive == 0
    upperbound, lowerbound = 0.1, 0.8
    bandwidth, numband = determine_band(Jaccard1, Jaccard2, upperbound, lowerbound)
    #bandwidth, numband = 5, 10
    @show upperbound, lowerbound
    lsh_vars = LSH_vars(numband, bandwidth)
    ret_prob_table(numband, bandwidth)
    # Following if
    println("sigmat_file not found generate sigmat...")
    sigmat = generate_signature_matrix(lsh_vars, binarray_csc, window_num=window_num, slidewidth=slidewidth)
    #@time sigmat = generate_signature_matrix(lsh_vars, binarray_csc, window_num=window_num, slidewidth=slidewidth)
    # open(io -> serialize(io, sigmat), sigmat_file, "w")
    # sigmat_file = joinpath(source_dir, "sigmat_mlenseq_"*string(min_len_sequence)*".dat")
    # if !ispath(sigmat_file)
    #     println("sigmat_file not found generate sigmat...")
    #     sigmat = generate_signature_matrix(lsh_vars, binarray_csc, window_num=window_num, slidewidth=slidewidth)
    #     #@time sigmat = generate_signature_matrix(lsh_vars, binarray_csc, window_num=window_num, slidewidth=slidewidth)
    #     open(io -> serialize(io, sigmat), sigmat_file, "w")
    # else
    #     println("sigmat_file found loading it...")
    #     sigmat = open(deserialize, sigmat_file)
    # end
    println("generate buckets_list")
    #@time buckets_list = generate_buckets(lsh_vars, sigmat);
    buckets_list = generate_buckets(lsh_vars, sigmat);
    candidates_set = Tuple{Int64, Int64}[]
    candidates = Set{Int64}()
    jaccard(set1, set2) = length(intersect(set1, set2)) / length(union(set1, set2))
    println("Start generation of candidates_set")
    @eval @everywhere times = $times
    @eval @everywhere sigmat = $sigmat
    @eval @everywhere lsh_vars = $lsh_vars
    @eval @everywhere buckets_list = $buckets_list
    @everywhere function pfind_similar(idx)
        idx1 = times[idx]
        indices = [times[idx2] for idx2 in find_similar(lsh_vars, sigmat, buckets_list, idx)]
        candidates_set = Tuple{Int64, Int64}[]    
        for idx2 in indices
            if idx1 != idx2
                push!(candidates_set, (idx1, idx2))
            end
        end
        return(candidates_set)
    end
    index_list = [idx for (idx, col) in enumerate(times)]
    println("start seeking for similar pairs")
    candidates_pair_ = pmap(pfind_similar, index_list);
    #@time candidates_pair_ = pmap(pfind_similar, index_list);
    # concatinate the result
    candidates_set = Tuple{Int64, Int64}[]
    for cp in candidates_pair_
        append!(candidates_set, cp)
    end
    candidates_set = unique(candidates_set);
    @eval @everywhere candidates_set = $candidates_set
    @everywhere jaccard(set1, set2) = length(intersect(set1, set2)) / length(union(set1, set2))
    @everywhere function find_most_similar_pair(idx1, idx2)
        m11 = Set{Int64}(findn(binarray_csc[:, idx1:(idx1+window_num-1)] .== 1)[1])
        m12 = Set{Int64}(findn(binarray_csc[:, (idx1+slidewidth):(idx1+slidewidth+window_num-1)] .== 1)[1])
        m21 = Set{Int64}(findn(binarray_csc[:, idx2:(idx2+window_num-1)] .== 1)[1])
        m22 = Set{Int64}(findn(binarray_csc[:, (idx2+slidewidth):(idx2+slidewidth+window_num-1)] .== 1)[1])
        # NeuronIDの一致が５つ以上ないようであればその計算は省く
        j1121 = jaccard(m11, m21)
        j1122 = jaccard(m11, m22)
        j1221 = jaccard(m12, m21)
        j1222 = jaccard(m12, m22)
        maxj = max(j1121, j1122, j1221, j1222)
        if (maxj == j1121)
            return(idx1, idx2)
        elseif (maxj == j1122)
            return(idx1, idx2+slidewidth)
        elseif (maxj == j1221)
            return(idx1+slidewidth, idx2)
        else
            return(idx1+slidewidth, idx2+slidewidth)
        end
    end
    candidates_set_alt = pmap(idx->find_most_similar_pair(idx[1], idx[2]), candidates_set)
    #@time candidates_set_alt = pmap(idx->find_most_similar_pair(idx[1], idx[2]), candidates_set)
else
    println("Exhaustive case")
    candidates_set_alt = Tuple{Int64, Int64}[]
    for t1 in times
        for t2 in times
            if  t1 != t2
                push!(candidates_set_alt, (t1, t2))
            end
        end
    end
    candidates_set_alt = unique(candidates_set_alt)
end
println("reduce rate", length(candidates_set_alt) / length(times)^2)
@eval @everywhere a = $a
@everywhere function pcal_editsim(row::Int, col::Int)
    data = edit_sim_local_with_gap(binarray_csc[:, row:(row+window_num-1)], binarray_csc[:, col:(col+window_num-1)], a=a, base=0)
    return([row, col, data])
end
@time results = pmap(idx->pcal_editsim(idx[1], idx[2]), candidates_set_alt)
#@time results = pmap(idx->pcal_editsim(idx[1], idx[2]), candidates_set_alt)
results_ = hcat(results...)';
@show results_
rows = results_[:, 1];
cols = results_[:, 2];
datas = results_[:, 3];
println("finished simmat calculation!! Now results are being stored...")
npzwrite(joinpath(save_dir, "parameters.npz"),
 Dict("a" => a,
      "window(ms)" => window,
      "slide(ms)" => slidewidth*binsize,
      "min_len_sequence" => min_len_sequence
     )
)
npzwrite(joinpath(save_dir, "simmat_coo.npz"),
  Dict("row" => rows - 1,
  "col" => cols - 1,
  "data" => datas,
  "dummy" => 1
  ))

println("finished simmat cal!!!")

