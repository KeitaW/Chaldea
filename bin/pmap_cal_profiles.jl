using NPZ
using PyCall
using MatrixMarket

@show nworkers(), nprocs()
@everywhere include("../chaldea/edit_sim.jl")
@everywhere include("../chaldea/barton_sternberg.jl")
@everywhere include("../chaldea/LSH.jl")
@everywhere include("../chaldea/utils/myutils.jl")
ret_dotdot = myutils.ret_dotdot
#function load_copra_data(rf; drop=3)
#    lines = open(rf, "r") do fp
#    readlines(fp)
#    end
#    data = []
#    for line in lines
#        tmp = [parse(Int, l) for l in split(line, ' ')[1:(end-1)]]
#        if length(tmp) >= drop
#            push!(data, tmp)
#        end
#    end
#    return(data)
#end

# originalSTDOUT = STDOUT
# (outRead, outWrite) = redirect_stdout()
script_dir = pwd()
source_dir = ARGS[1]
numiter = parse(Int64, ARGS[2]) # unit less
min_num_seq = parse(Int64, ARGS[3])

@show source_dir
@show numiter
rf = joinpath(source_dir, "best-clusters-rereduced_simmat.mtx")
data = myutils.load_copra_data(rf)
@show npzread(joinpath(source_dir, "correspondence_table_time.npz"))
correspondence_table_time = npzread(joinpath(source_dir, "correspondence_table_time.npz"))["arr_0"]
reference_points_list = []
for d in data
    push!(reference_points_list, correspondence_table_time[d] + 1)
end
# clusters dir
parameters = npzread(joinpath(ret_dotdot(source_dir), "parameters.npz"))

a = parameters["a"]
window = parameters["window(ms)"]
slidewidth = parameters["slide(ms)"]
d = chomp(readstring(`date '+%Y%m%dT%H%M%S'`))
save_dir = joinpath(source_dir, "profiles_"*"numiter_"*string(numiter)*"_"*d)
@show save_dir
if !isdir(save_dir)
    mkdir(save_dir)
end
open(fid->write(fid, "$save_dir"), "/tmp/save_dir", "w")
binarray_data = npzread(joinpath(ret_dotdot(source_dir, 2), "binarray_data.npz"))
rows = binarray_data["row"] + 1
cols = binarray_data["col"] + 1
datas = binarray_data["data"]
binsize = binarray_data["binsize(ms)"]
duration = binarray_data["duration(ms)"]
window_num = div(window, binsize)
@show window, window_num
rows = [Int64(r) for r in rows];
cols = [Int64(r) for r in cols];
binarray_csc = sparse(rows, cols, datas)
@show binsize, size(binarray_csc)
# reference_pointsに対応する部分神経活動の抽出
@everywhere mats_list = []
for reference_points in reference_points_list
    mats = []
    for (idx, reference_timing) in enumerate(reference_points)
        rt = Int(reference_timing)
        if rt >= size(binarray_csc)[2] - window_num
            continue
        end
        # 最低限3回の発火が無いと解析では扱わない．のけておかないとProfileがNullってしまう
        if sum(binarray_csc[:, rt:(rt+window_num)]) <= 3
            continue
        else
            @show size(full(binarray_csc[:, rt:(rt+window_num)]))
            push!(mats, full(binarray_csc[:, rt:(rt+window_num)]))
            @show length(mats)
        end
    end
    if length(mats) >= 1
        push!(mats_list, mats)
    end
end
@eval @everywhere numiter = $numiter
@eval @everywhere window_num = $window_num
@eval @everywhere a = $a
@everywhere function pret_profile(idx, mats)
    if length(mats) >= 2
        bs_pre = BartonSternberg(mats, window_num, maxiter = numiter, a = a, base = 0)
        bs_post = align_mats(bs_pre)
        return(Dict{String, Any}(string(idx) => bs_post.profile))
    else
        return nothing
    end

end

for (idx, mats) in enumerate(mats_list)
    @show pret_profile(idx, mats)
    for mat in mats
        @show size(mat)
    end
end
profiles_dict = merge(pmap(x->pret_profile(x[1], x[2]), enumerate(mats_list))...)
@show profiles_dict
npzwrite(joinpath(save_dir, "profiles.npz"), profiles_dict)
