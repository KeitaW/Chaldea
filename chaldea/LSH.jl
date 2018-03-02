using  ProgressMeter
# こちらのソースファイルの内容を変更する際には../Notes/12_17_test_MINLSH.ipynbのコードを変更・テストの後コピペすること
type LSH_vars
    # LSHに用いる各種変数を格納しておくtype
    numhash::Int64
    numband::Int64
    bandwidth::Int64
    function LSH_vars(nb, bw)
        new(nb*bw, nb, bw)
    end
end

# 多重ディスパッチ
function minhash(words::Array{Int64, 1}, seed::UInt64)
    current_min = typemax(UInt64)
    minhash_word = nothing
    for word in words
        hash_ = hash(word, seed)
        if hash_ < current_min
            minhash_word = word
            current_min = hash_
        end
    end
    return(minhash_word)
end
function minhash(words::Set{Int64}, seed::UInt64)
    current_min = typemax(UInt64)
    minhash_word = nothing
    for word in words
        hash_ = hash(word, seed)
        if hash_ < current_min
            minhash_word = word
            current_min = hash_
        end
    end
    return(minhash_word)
end

function generate_signature_matrix(lsh_vars::LSH_vars, binarray_csc::SparseMatrixCSC{Float64,Int64}; window_num=100, slidewidth=50)
    numhash = lsh_vars.numhash
    extent = 1:slidewidth:(size(binarray_csc)[2] - 2window_num)
    signature_matrix = zeros((numhash, length(extent)))
    p = Progress(numhash, 1)
    for row in 1:numhash
        for (idx, col) in enumerate(extent)
            # 発火が存在するIDのみからなるSetを作成する
            idsets = Set{Int64}(findn(binarray_csc[:, col:(col+window_num-1)] .== 1)[1])
            if length(idsets) > 0
                signature_matrix[row, idx] = minhash(idsets, UInt64(row))
            else
                # idsetsが空の場合 = 発火が存在しない場合も衝突を回避するためランダムな値を詰めておく
                signature_matrix[row, idx] = minhash(rand(1:typemax(Int64), 1), UInt64(row))
            end
         end
         next!(p)
    end
    return(signature_matrix)
end

# b: numband, r:bandwidth(rows)
cal_prob(s, b, r) = 1 - ( 1- s^r )^b
function ret_prob_table(b, r, step=0.1)
    for s in 0:step:1
        println(s, ' ', cal_prob(s, b, r))
    end
end

function generate_buckets(lsh_vars::LSH_vars, signature_matrix::Array{Float64,2})
    # signature matrixの各列をバケツに放り込む関数
    # 似ているものは同じバケツに放り込まれやすい．具体的な確率はバンドの設計による
    numhash = lsh_vars.numhash
    numband = lsh_vars.numband
    bandwidth = lsh_vars.bandwidth
    buckets_list = []
    for band in 1:bandwidth:numhash
        buckets = Dict{UInt64, Set{Int64}}()
        for col in 1:size(signature_matrix)[2]
            hash_ = hash(signature_matrix[band:(band+bandwidth-1), col], UInt64(band))
            if !haskey(buckets, hash_)
                buckets[hash_] = Set{Int64}()
            end
            push!(buckets[hash_], col)
        end
        push!(buckets_list, buckets)
    end
    return(buckets_list)
end

function find_similar(lsh_vars::LSH_vars, signature_matrix::Array{Float64,2}, buckets_list::Array{Any,1}, col::Int64)
    # buckets_listの情報を元にcolに類似した（≃ Jaccard係数が高い）要素を抽出するための関数
    # Python版とは異なり，編集類似度を用いたスクリーニングはここでは行わない
    numhash = lsh_vars.numhash
    numband = lsh_vars.numband
    bandwidth = lsh_vars.bandwidth
    candidates = Set{Int64}()
    for (idx, band) in enumerate(1:numhash:bandwidth)
        hash_ = hash(signature_matrix[band:(band+bandwidth-1), col], UInt64(band))
        candidates = union(candidates, buckets_list[idx][hash_])
    end
    return(candidates)
end

