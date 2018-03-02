# こちらのソースファイルの内容を変更する際には../Notes/12_16_test_local_alignment_with_gap.ipynbのコードを変更・テストの後コピペすること
# 一応σとbaseは分けたがこれらは同一でもいいかもしれない．
# また，ρはジッター幅を考えてbaseの値がその幅を超えると0になるように調整する．
function edit_sim_local_with_gap(m1::Array{Float64, 2}, m2::Array{Float64, 2}; base=0.01, a=0.1)
    # ±jitter を許す．ただし，線形でギャップペナルティをかける．
    # ギャップを導入すること自体にはペナルティをかけていないことに注意
    # ρはギャップを導入するペナルティ，σはギャップを伸ばすペナルティ
    # baseは内積計算の際のベースラインを決定するもの．最初にbinarrayから列をサンプリングして決める
    # aはギャップの伸ばしていった際に指数関数的にペナルティを欠ける際の重み項
    mat1 = flipdim(m1, 2)
    mat2 = flipdim(m2, 2)
    ncol1 = size(mat1)[2]
    ncol2 = size(mat2)[2]
    dp = zeros(Float64, (ncol1+1, ncol2+1))
    # affine gapの要領でいくつギャップを伸ばしたかを記録する
    down = zeros(Int64, (ncol1+1, ncol2+1))
    right = zeros(Int64, (ncol1+1, ncol2+1))
    for col1 in 1:ncol1
        for col2 in 1:ncol2
            match = dot(mat1[:, col1],  mat2[:, col2])
            # downの判定 新しくギャップを伸ばす vs 更にギャップを伸ばす
            which_dir = indmax([dp[col1, col2+1] - exp(a)+exp(0),
                dp[col1-down[col1, col2+1], col2+1] - exp(a*(down[col1, col2+1]+1)) + exp(0)])
            if which_dir == 1
                down[col1+1, col2+1] = 1
            else
                down[col1+1, col2+1] = down[col1, col2+1] + 1
            end
            down_score = dp[col1-down[col1+1, col2+1]+1, col2+1] - exp(a*down[col1+1, col2+1]) + exp(0)
            # rightの判定 新しくギャップを伸ばす vs 更にギャップを伸ばす
            which_dir = indmax([dp[col1+1, col2] - exp(a)  + exp(0),
                dp[col1+1, col2-right[col1+1, col2]] - exp(a*right[col1+1, col2]+1) + exp(0)])
            if which_dir == 1
                right[col1+1, col2+1] = 1
            else
                right[col1+1, col2+1] = right[col1+1, col2] + 1
            end
            right_score = dp[col1+1, col2-right[col1+1, col2+1]+1] - exp(a*right[col1+1, col2+1]) + exp(0)
            # 最初の0は始点から任意の頂点への"タダ乗り"．ローカル配列アラインメントを求めている
            dp[col1+1, col2+1] = max(0, down_score , right_score, dp[col1, col2] + match - base)
        end
    end
    # 実質的には任意の点から終点への"タダ乗り"
    return(maximum(dp))
end
function edit_sim_local_with_gap(m1::SparseMatrixCSC{Float64,Int64},
    m2::SparseMatrixCSC{Float64,Int64}; base=0.01, a=0.1)
    # ±jitter を許す．ただし，線形でギャップペナルティをかける．
    # ギャップを導入すること自体にはペナルティをかけていないことに注意
    # ρはギャップを導入するペナルティ，σはギャップを伸ばすペナルティ
    # baseは内積計算の際のベースラインを決定するもの．最初にbinarrayから列をサンプリングして決める
    # aはギャップの伸ばしていった際に指数関数的にペナルティを欠ける際の重み項
    mat1 = flipdim(m1, 2)
    mat2 = flipdim(m2, 2)
    ncol1 = size(mat1)[2]
    ncol2 = size(mat2)[2]
    dp = zeros(Float64, (ncol1+1, ncol2+1))
    # affine gapの要領でいくつギャップを伸ばしたかを記録する
    down = zeros(Int64, (ncol1+1, ncol2+1))
    right = zeros(Int64, (ncol1+1, ncol2+1))
    for col1 in 1:ncol1
        for col2 in 1:ncol2
            match = sum(mat1[:, col1] .* mat2[:, col2])
            # downの判定 新しくギャップを伸ばす vs 更にギャップを伸ばす
            which_dir = indmax([dp[col1, col2+1] - exp(a) + exp(0),
                dp[col1-down[col1, col2+1], col2+1] - exp(a*(down[col1, col2+1]+1)) + exp(0)])
            if which_dir == 1
                down[col1+1, col2+1] = 1
            else
                down[col1+1, col2+1] = down[col1, col2+1] + 1
            end
            down_score = dp[col1-down[col1+1, col2+1]+1, col2+1] - exp(a*down[col1+1, col2+1]) + exp(0)
            # rightの判定 新しくギャップを伸ばす vs 更にギャップを伸ばす
            which_dir = indmax([dp[col1+1, col2] - exp(a),
                dp[col1+1, col2-right[col1+1, col2]] - exp(a*right[col1+1, col2]+1)]) + exp(0)
            if which_dir == 1
                right[col1+1, col2+1] = 1
            else
                right[col1+1, col2+1] = right[col1+1, col2] + 1
            end
            right_score = dp[col1+1, col2-right[col1+1, col2+1]+1] - exp(a*right[col1+1, col2+1]) + exp(0)
            # 最初の0は始点から任意の頂点への"タダ乗り"．ローカル配列アラインメントを求めている
            dp[col1+1, col2+1] = max(0, down_score , right_score, dp[col1, col2] + match - base)
        end
    end
    # 実質的には任意の点から終点への"タダ乗り"
    return(maximum(dp))
end

function edit_sim_local_with_gap_with_bp(m1::Array{Float64, 2}, m2::Array{Float64, 2}; ρ=1, σ=1, base=0.01, a=0.05)
  # ±jitter を許す．ただし，線形でギャップペナルティをかける．
  # ギャップを導入すること自体にはペナルティをかけていないことに注意
  # ρはギャップを導入するペナルティ，σはギャップを伸ばすペナルティ
  # baseは内積計算の際のベースラインを決定するもの．最初にbinarrayから列をサンプリングして決める
  mat1 = flipdim(m1, 2)
  mat2 = flipdim(m2, 2)
  ncol1 = size(mat1)[2]
  ncol2 = size(mat2)[2]
  dp = zeros(Float64, (ncol1+1, ncol2+1))
  bp = ones(Int64, (ncol1+1, ncol2+1)) * (-1)
  down = zeros(Int64, (ncol1+1, ncol2+1))
  right = zeros(Int64, (ncol1+1, ncol2+1))
  for col1 in 1:ncol1
      for col2 in 1:ncol2
          match = dot(mat1[:, col1],  mat2[:, col2])
          # downの判定 新しくギャップを伸ばす vs 更にギャップを伸ばす
          which_dir = indmax([dp[col1, col2+1] - exp(a) + exp(0),
              dp[col1-down[col1, col2+1], col2+1] - exp(a*(down[col1, col2+1]+1)) + exp(0)])
          if which_dir == 1
              down[col1+1, col2+1] = 1
          else
              down[col1+1, col2+1] = down[col1, col2+1] + 1
          end
          down_score = dp[col1-down[col1+1, col2+1]+1, col2+1] - exp(a*down[col1+1, col2+1]) + exp(0)
          # rightの判定 新しくギャップを伸ばす vs 更にギャップを伸ばす
          which_dir = indmax([dp[col1+1, col2] - exp(a) + exp(0),
              dp[col1+1, col2-right[col1+1, col2]] - exp(a*right[col1+1, col2]+1) + exp(0)])
          if which_dir == 1
              right[col1+1, col2+1] = 1
          else
              right[col1+1, col2+1] = right[col1+1, col2] + 1
          end
          right_score = dp[col1+1, col2-right[col1+1, col2+1]+1] - exp(a*right[col1+1, col2+1]) + exp(0)
          # 最初の0は始点から任意の頂点への"タダ乗り"．ローカル配列アラインメントを求めている
          tmp = [0, down_score, right_score, dp[col1, col2] + match - base]
          maxdir = indmax(tmp)
          bp[col1+1, col2+1] = maxdir
          dp[col1+1, col2+1] = tmp[maxdir]
      end
  end
  # 通常であればここはdp[end, end]でよい．ただし今回のようなギャップペナルティを導入すると
  # 共通配列が各配列の最初の方にあった場合，そこでプラスのスコアを記録しても後のギャップ導入で
  # 結局スコアが0ないしマイナスになってしまう問題が考えられる．
  # そこでその問題を解決するべくdpの最大値を返すように変更を加える
  # さらに，後の共通配列計算のため，dpの最大値に対応する場所を返してやるようにする
  # このmaxrow, maxcolの位置からalignmentを始めないとtragedy
  maxrow, maxcol = ind2sub(size(dp), indmax(dp))
  return(bp, dp, maxrow, maxcol)
end
function edit_sim_local_with_gap_with_bp(m1::SparseMatrixCSC{Float64,Int64},
    m2::SparseMatrixCSC{Float64,Int64}; base=0.01, a=0.1)
    # ±jitter を許す．ただし，線形でギャップペナルティをかける．
    # ギャップを導入すること自体にはペナルティをかけていないことに注意
    # ρはギャップを導入するペナルティ，σはギャップを伸ばすペナルティ
    # baseは内積計算の際のベースラインを決定するもの．最初にbinarrayから列をサンプリングして決める
    mat1 = flipdim(m1, 2)
    mat2 = flipdim(m2, 2)
    ncol1 = size(mat1)[2]
    ncol2 = size(mat2)[2]
    dp = zeros(Float64, (ncol1+1, ncol2+1))
    bp = ones(Int64, (ncol1+1, ncol2+1)) * (-1)
    down = zeros(Int64, (ncol1+1, ncol2+1))
    right = zeros(Int64, (ncol1+1, ncol2+1))
    for col1 in 1:ncol1
        for col2 in 1:ncol2
            match = sum(mat1[:, col1] .* mat2[:, col2])
            # downの判定 新しくギャップを伸ばす vs 更にギャップを伸ばす
            which_dir = indmax([dp[col1, col2+1] - exp(a) + exp(0),
                dp[col1-down[col1, col2+1], col2+1] - exp(a*(down[col1, col2+1]+1)) + exp(0)])
            if which_dir == 1
                down[col1+1, col2+1] = 1
            else
                down[col1+1, col2+1] = down[col1, col2+1] + 1
            end
            down_score = dp[col1-down[col1+1, col2+1]+1, col2+1] - exp(a*down[col1+1, col2+1]) + exp(0)
            # rightの判定 新しくギャップを伸ばす vs 更にギャップを伸ばす
            which_dir = indmax([dp[col1+1, col2] - exp(a) + exp(0),
                dp[col1+1, col2-right[col1+1, col2]] - exp(a*right[col1+1, col2]+1) + exp(0)])
            if which_dir == 1
                right[col1+1, col2+1] = 1
            else
                right[col1+1, col2+1] = right[col1+1, col2] + 1
            end
            right_score = dp[col1+1, col2-right[col1+1, col2+1]+1] - exp(a*right[col1+1, col2+1]) + exp(0)
            # 最初の0は始点から任意の頂点への"タダ乗り"．ローカル配列アラインメントを求めている
            tmp = [0, down_score, right_score, dp[col1, col2] + match - base]
            maxdir = indmax(tmp)
            bp[col1+1, col2+1] = maxdir
            dp[col1+1, col2+1] = tmp[maxdir]
        end
    end
    # 通常であればここはdp[end, end]でよい．ただし今回のようなギャップペナルティを導入すると
    # 共通配列が各配列の最初の方にあった場合，そこでプラスのスコアを記録しても後のギャップ導入で
    # 結局スコアが0ないしマイナスになってしまう問題が考えられる．
    # そこでその問題を解決するべくdpの最大値を返すように変更を加える
    # さらに，後の共通配列計算のため，dpの最大値に対応する場所を返してやるようにする
    # このmaxrow, maxcolの位置からalignmentを始めないとtragedy
    maxrow, maxcol = ind2sub(size(dp), indmax(dp))
    return(bp, dp, maxrow, maxcol)
end
# bp tableの情報を元に2つの神経活動をalignする関数
function align_editsim(bp, maxrow, maxcol, m1, m2)
    mat1 = flipdim(m1, 2)
    mat2 = flipdim(m2, 2)
    zerovec = zeros(Float64, (size(mat1)[1], 1))
    # 先頭は単純に初期化エラーを回避するために入れている．返す際には無視すること
    alignment1 = zeros(Float64, (size(mat1)[1], 1))
    alignment2 = zeros(Float64, (size(mat1)[1], 1))
    i = size(bp)[1]
    j = size(bp)[2]
   # 終端へのジャンプも導入しているため，まずdpの最大値に対応する地点までi, jを動かす
#    while i > maxrow
#    alignment1 = cat(2, mat1[:, i-1], alignment1)
#    i -= 1
#    end
#    while j > maxcol
#        alignment2 = cat(2, mat2[:, j-1], alignment2)
#        j -= 1
#    end
    # ジャンプを無視する場合
    #i = maxrow
    #j = maxcol
    while true
        if i == 1 || j == 1
            # for test!!!
            break
        elseif bp[i, j] == 4
            # 一致があったことに対応
            if sum(mat1[:, i-1]) != 0 && sum(mat2[:, j-1] != 0)
                alignment1 = cat(2, round(mat1[:, i-1] .* mat2[:, j-1]), alignment1)
                alignment2 = cat(2, round(mat1[:, i-1] .* mat2[:, j-1]), alignment2)
            end
            i -= 1
            j -= 1
        # # 上の別バージョン
        # elseif bp[i, j] == 4
        #     # 一致があったことに対応
        #     if sum(mat1[:, i-1]) != 0 && sum(mat2[:, j-1] != 0)
        #         alignment1 = cat(2, mat1[:, i-1], alignment1)
        #         alignment2 = cat(2, mat2[:, j-1], alignment2)
        #     end
        #     i -= 1
        #     j -= 1
        elseif bp[i, j] == 2
            # down
            alignment1 = cat(2, mat1[:, i-1], alignment1)
            alignment2 = cat(2, zerovec, alignment2)
            i -= 1
        elseif bp[i, j] == 3
            # right
            alignment1 = cat(2, zerovec, alignment1)
            alignment2 = cat(2, mat2[:, j-1], alignment2)
            j -= 1
        elseif bp[i, j] == -1
            # 端に至ったことに対応する
            break
        elseif bp[i, j] == 1
            # 始点への高飛びを表す
            # return (flipdim(alignment1, 2), flipdim(alignment2, 2))
            # そのまま返してしまうと最終的に行列が消滅してしまうため，以下のコードをテスト(2015/12/27 11:52)
            break
        end
    end
    # 始点への高飛び後の後処理
    if i >= 2
        for remaini in i:-1:2
            alignment1 = cat(2, mat1[:, remaini-1], alignment1)
        end
    end
    if j >= 2
        for remainj in j:-1:2
            alignment2 = cat(2, mat2[:, remainj-1], alignment2)
        end
    end
    # cat関数を用いているためか，sparse -> floatの変換が行われてしまっている点に注意.
    # 1:end-1: 初期の分を除いている
    return (flipdim(alignment1[:, 1:end-1], 2), flipdim(alignment2[:, 1:end-1], 2))
end

# bp tableの情報を元に共通成分を発見する関数m1が神経活動,m2がprofileとなる

# bp tableの情報を元に2つの神経活動をalignする関数
function find_common(bp, maxrow, maxcol, m1, m2; hosei=0.3)
    mat1 = flipdim(m1, 2)
    mat2 = flipdim(m2, 2)
    zerovec = zeros(Float64, (size(mat1)[1], 1))
    # 先頭は単純に初期化エラーを回避するために入れている．返す際には無視すること
    alignment1 = zeros(Float64, (size(mat1)[1], 1))
    alignment2 = zeros(Float64, (size(mat1)[1], 1))
    i = size(bp)[1]
    j = size(bp)[2]
   # 終端へのジャンプも導入しているため，まずdpの最大値に対応する地点までi, jを動かす
    while i > maxrow
    alignment1 = cat(2, zerovec, alignment1)
    i -= 1
    end
    while j > maxcol
        alignment2 = cat(2, zerovec, alignment2)
        j -= 1
    end
    # ジャンプを無視する場合
    #i = maxrow
    #j = maxcol
    while true
        if i == 1 || j == 1
            # for test!!!
            break
        elseif bp[i, j] == 4
            # 一致があったことに対応
            alignment1 = cat(2, round(mat1[:, i-1] .* mat2[:, j-1] .+ hosei), alignment1)
            alignment2 = cat(2, round(mat1[:, i-1] .* mat2[:, j-1] .+ hosei), alignment2)
            i -= 1
            j -= 1
        elseif bp[i, j] == 2
            # down
            alignment1 = cat(2, zerovec, alignment1)
            # 変更
            # alignment1 = cat(2, mat1[:, i-1], alignment1)
            # alignment2 = cat(2, zerovec, alignment2)
            i -= 1
        elseif bp[i, j] == 3
            # right
            # alignment1 = cat(2, zerovec, alignment1)
            alignment2 = cat(2, zerovec, alignment2)
            j -= 1
        elseif bp[i, j] == -1
            # 端に至ったことに対応する
            break
        elseif bp[i, j] == 1
            # 始点への高飛びを表す
            # return (flipdim(alignment1, 2), flipdim(alignment2, 2))
            # そのまま返してしまうと最終的に行列が消滅してしまうため，以下のコードをテスト(2015/12/27 11:52)
            break
        end
    end
    # 始点への高飛び後の後処理
    if i >= 2
        for remaini in i:-1:2
            alignment1 = cat(2, zerovec, alignment1)
        end
    end
    if j >= 2
        for remainj in j:-1:2
            alignment2 = cat(2, zerovec, alignment2)
        end
    end
    # cat関数を用いているためか，sparse -> floatの変換が行われてしまっている点に注意.
    # 1:end-1: 初期の分を除いている
    return (flipdim(alignment1[:, 1:end-1], 2), flipdim(alignment2[:, 1:end-1], 2))
end
