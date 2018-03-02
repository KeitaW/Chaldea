using Distributions
using ProgressMeter
using NPZ
include("edit_sim.jl")
type BartonSternberg
    mats::Array{Array{Float64,2},1}
    # 元の行列を保持する．毎回これをalignmentしてやること．matsはcolが欠落するため
    orgmats::Array{Array{Float64,2},1}
    profile::Array{Float64,2}
    nmats::Int64
    unprosessed::Array{Bool, 1}
    processed_ids::Array{Int64, 1}
    numiter::Int64
    maxiter::Int64
    a::Float64
    base::Float64
    numwindow::Int64
    function BartonSternberg(mats, numwindow; maxiter=10, a=0.05, base=0.01)
        profile = zeros(Float64, size(mats[1]))
        new(mats[:], mats[:], profile, length(mats), ones(Bool, length(mats)), zeros(Int64, length(mats)), 0, maxiter, 0, 0, numwindow)
    end
end


function apply_gaussian_filter(mat, σ)
    if σ == 0
        return mat
    end
    @show size(mat)
    normal = Normal(0, σ)
    extent = 3σ
    normal_filter = transpose(pdf(normal, -float(extent):float(extent)))
    rows, cols = findn(mat .> 0)
    filtered_mat = zeros(size(mat)[1], (size(mat)[2]+2extent))
    for (row, col) in zip(rows, cols)
      # はみ出し防止のため，ちょっと大きめに行列を確保している
      # colのズレを修正
      col += extent
      @show size(filtered_mat[row, col-extent:col+extent]), size(normal_filter)
      filtered_mat[row, col-extent:col+extent] += normal_filter[1:2extent]
    end
    # "みみ"の部分を切り捨てる
    filtered_mat = filtered_mat[:, (extent+1):(end-extent)]
    return(filtered_mat)
end

function shrink_profile(profile, numwindow)
    # col to col to width
    ncol = size(profile)[2]
    if ncol <= numwindow
        return(profile)
    end
    ctoc_width = round(Int, ncol / numwindow)
    ctoc_width = ctoc_width % 2 == 0 ? ctoc_width : ctoc_width + 1
    @show ncol, numwindow, ctoc_width
    span = div(ctoc_width, 2):ctoc_width:size(profile)[2]-2div(ctoc_width, 2)
    shrinked_profile = zeros(size(profile)[1], length(span));
    for row in 1:size(profile)[1]
        for (shrinked_col, col) in enumerate(span)
            col += 1
            shrinked_profile[row, shrinked_col] = sum(profile[row, (col-div(ctoc_width, 2)):(col+div(ctoc_width, 2))]) / ctoc_width
        end
    end
    return(shrinked_profile)
end

function eliminate_empty_colums(mat)
    # すべての要素が0の列を省く処理
    indices_remain = Int64[]
    rowsums = sum(mat, 1)
    for col in 1:size(mat)[2]
        if rowsums[col] != 0
            push!(indices_remain, col)
        end
    end
    mat = mat[:, indices_remain]
    return(mat)
end

function regularize_row(mat)
    for row in 1:size(mat)[1]
        if (maximum(mat[row, :]) - minimum(mat[row, :])) != 0
            mat[row, :] = (mat[row, :] .- mean(mat[row, :])) ./ var(mat[row, :])
            mat[row, :] = (mat[row, :] .- minimum(mat[row, :])) ./ (maximum(mat[row, :]) .- minimum(mat[row, :]))
        end
    end
    return(mat)
end

function _generate_profile(bs::BartonSternberg, σ::Int64)
    nrow = size(bs.mats[1])[1]
    ncol = maximum([size(mat)[2] for mat in bs.mats])
    cummat = zeros(nrow, ncol)
    for mat in bs.mats
        cummat[:, 1:size(mat)[2]] += mat[:, :]
    end
    cummat = eliminate_empty_colums(cummat)
    cummat = regularize_row(cummat)
    # bs.profile = apply_gaussian_filter(cummat, σ)
    bs.profile = cummat
end

# 多重ディスパッチ．idxに対応する行列を除いた形でprofileの作成を行う
function _generate_profile(bs::BartonSternberg, idx::Int64, σ::Int64)
    # リスト内包表記で条件式が使えないため，こういう形で書いている
    indices = ones(Bool, bs.nmats)
    indices[idx] = false
    mats = bs.mats[indices]
    nrow = size(mats[1])[1]
    ncol = maximum([size(mat)[2] for mat in mats])
    cummat = zeros(nrow, ncol)
    for mat in mats
        cummat[:, 1:size(mat)[2]] += mat[:, :]
    end
    cummat = eliminate_empty_colums(cummat)
    cummat = regularize_row(cummat)
    # bs.profile = apply_gaussian_filter(cummat, σ)
    bs.profile = cummat
end

# 多重ディスパッチ．indicesに対応する行列のみでprofileの作成を行う
function _generate_profile(bs::BartonSternberg, indices::Array{Bool, 1}, σ::Int64)
    # リスト内包表記で条件式が使えないため，こういう形で書いている
    mats = bs.mats[indices]
    nrow = size(mats[1])[1]
    ncol = maximum([size(mat)[2] for mat in mats])
    cummat = zeros(nrow, ncol)
    for mat in mats
        cummat[:, 1:size(mat)[2]] += mat[:, :]
    end
    cummat = eliminate_empty_colums(cummat)
    cummat = regularize_row(cummat)
    # bs.profile = apply_gaussian_filter(cummat, σ)
    bs.profile = cummat
end


# 多重ディスパッチ．indicesに対応する行列のみでprofileの作成を行う
function _generate_profile(bs::BartonSternberg, indices::BitArray{1}, σ::Int64)
    # リスト内包表記で条件式が使えないため，こういう形で書いている
    mats = bs.mats[indices]
    nrow = size(mats[1])[1]
    ncol = maximum([size(mat)[2] for mat in mats])
    cummat = zeros(nrow, ncol)
    for mat in mats
        cummat[:, 1:size(mat)[2]] += mat[:, :]
    end
    cummat = eliminate_empty_colums(cummat)
    cummat = regularize_row(cummat)
    # bs.profile = apply_gaussian_filter(cummat, σ)
    bs.profile = cummat
end

function _first_profile_generation(bs::BartonSternberg, σ::Int64)
#    # 時間がかかりすぎる
#    simmat = zeros(Float64, bs.nmats, bs.nmats)
#    for i in 1:(bs.nmats)
#       for j in 1:(i-1)
#           simmat[i, j] = edit_sim_local_with_gap(bs.orgmats[i], bs.orgmats[j], a=bs.a)
#       end
#    end
#    # i, j: most similar pairs
#    # まずはこの２つの添字に対応した行列をalginする
#    i, j = ind2sub(size(simmat), indmax(simmat))
    # # 上が時間を食い過ぎる場合の代替案
    i, j = 1, 2
    bs.unprosessed[i] = false
    bs.unprosessed[j] = false
    bs.processed_ids[i] =  bs.numiter
    bs.numiter += 1
    bs.processed_ids[j] =  bs.numiter
    bs.numiter += 1
    @show size(bs.orgmats[i]),size(bs.orgmats[j])
    bp, db, maxrow, maxcol = edit_sim_local_with_gap_with_bp(bs.orgmats[i], bs.orgmats[j], a=bs.a, base=bs.base)
    bs.mats[i], bs.mats[j] = align_editsim(bp, maxrow, maxcol, bs.orgmats[i], bs.orgmats[j])
    @show bs.unprosessed, typeof(bs.unprosessed)
    @show .!(bs.unprosessed), typeof(.!(bs.unprosessed))
    _generate_profile(bs, .!(bs.unprosessed), σ)
end

function _iterative_update(bs::BartonSternberg, σ::Int64, a::Float64)
  if sum(bs.unprosessed) != 0
        indices = collect(1:bs.nmats)[bs.unprosessed]
        maxsim = -1.0
        maxidx = 1
        for (idx, mat) in zip(indices, bs.orgmats[bs.unprosessed])
            sim = edit_sim_local_with_gap(mat, bs.profile, a=bs.a)
            if maxsim < sim
                # profileと最も近い行列の添字を格納する
                maxsim = sim
                maxidx = idx
            end
        end
        idx = maxidx
        bs.unprosessed[idx] = false
        bs.processed_ids[idx] = bs.numiter
        bs.numiter += 1
        # パターン2 profileが先に来る場合
        bp, db, maxrow, maxcol = edit_sim_local_with_gap_with_bp(bs.profile, bs.orgmats[idx], a=a, base=bs.base)
        _, bs.mats[idx] = align_editsim(bp, maxrow, maxcol, bs.profile, bs.orgmats[idx])
        # # パターン１ matが先に来る場合
        # bp, db, maxrow, maxcol = edit_sim_local_with_gap_with_bp(bs.orgmats[idx], bs.profile, a=a, base=bs.base)
        # bs.mats[idx], _ = align_editsim(bp, maxrow, maxcol, bs.orgmats[idx], bs.profile)
        # bp, db = edit_sim_local_with_gap_with_bp(bs.profile, bs.mats[idx], ρ=bs.ρ, σ=bs.σ)
        # _, bs.mats[idx] = align_editsim(bp, bs.profile, bs.mats[idx])
        # 既に処理が終了したものだけでprofileを作成する
        _generate_profile(bs, !bs.unprosessed, σ)
    else
        # すべての行列をなめ終わったら，同様の順番でprofileを更新していく
        idx = indmin(bs.processed_ids)
        bs.numiter += 1
        bs.processed_ids[idx] = bs.numiter

        # # (上の代替案) すべての行列をなめ終わったら，ランダムな順番でprofileを更新していく
        # idx = sample(1:length(bs.mats))
        # bs.numiter += 1
        # bs.processed_ids[idx] = bs.numiter

        # idxの行列だけ除いてprofileを作成し，それを元に更新する
        _generate_profile(bs, idx)
        # パターン2 profileが先に来る場合
        bp, db, maxrow, maxcol = edit_sim_local_with_gap_with_bp(bs.profile, bs.orgmats[idx], a=a, base=bs.base)
        _, bs.mats[idx] = align_editsim(bp, maxrow, maxcol, bs.profile, bs.orgmats[idx])
        # # パターン１ matが先に来る場合
        # bp, db, maxrow, maxcol = edit_sim_local_with_gap_with_bp(bs.orgmats[idx], bs.profile, a=a, base=bs.base)
        # bs.mats[idx], _ = align_editsim(bp, maxrow, maxcol, bs.orgmats[idx], bs.profile)
        # bp, db = edit_sim_local_with_gap_with_bp(bs.profile, bs.mats[idx], ρ=bs.ρ, σ=bs.σ)
        # _, bs.mats[idx] = align_editsim(bp, bs.profile, bs.mats[idx])
        _generate_profile(bs, σ)
    end
end

function align_mats(bs::BartonSternberg; save_dir = "", idx = 1)
    # σ: フィルタの分散 最小値1まで持ってく
    # σ = div(bs.maxiter, 10)
    σ = 10
    max_sigma = σ
    _first_profile_generation(bs, σ)
    # progress = Progress(bs.maxiter)
    sigma_step = max_sigma / bs.maxiter
    for i in 1:bs.maxiter
        # σ = div(bs.maxiter - i, 10) + 1 # sigma == 0を避けるため
        σ = Int(floor(max_sigma - i*sigma_step))
        # @show a, σ
        _iterative_update(bs, σ, bs.a)
        # 1000回に１回，現時点での処理結果を辞書に格納
#        @show i
#        if i % 100 == 0
#            save_data = Dict{ASCIIString,Any}()
#            save_data["profile"*string(idx)] = bs.profile
#            npzwrite(joinpath(save_dir, "profiles_idx"*string(idx)*"_niter"*string(i)*".npz"), save_data)
#        end
        # bs.profile = shrink_profile(bs.profile, bs.numwindow)
        # next!(progress)
    end
    # bs.profile = (bs.profile .- minimum(bs.profile)) ./ (maximum(bs.profile) - minimum(bs.profile))
    return bs
end
