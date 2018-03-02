#-*-coding:utf-8-*-

import numpy as np
import sys
import matplotlib.pyplot as plt
from ccal_sim import *

def set_trace():
    from IPython.core.debugger import Pdb
    Pdb(color_scheme='Linux').set_trace(sys._getframe().f_back)

def debug(f, *args, **kwargs):
    from IPython.core.debugger import Pdb
    pdb = Pdb(color_scheme='Linux')
    return pdb.runcall(f, *args, **kwargs)


class BartonSternberg(object):
    """
    Barton-Sternberg multiple alignmentに際して，複数の行列を管理するためのクラス

    本アルゴリズムでalignを行った際には行列サイズが変わり得ることに注意
    """

    def __init__(self, mats, score_vec, jitter=20, window=400, max_iter=1000):
        self.mats           = list(mats)
        self.nmats          = len(mats)
        self.processed      = np.zeros(self.nmats, dtype=bool)
        self.processed_ids  = np.zeros(self.nmats, dtype=np.int)
        self.processed_id   = 0
        self.profile        = None
        self.lcss           = []
        self.max_iter       = max_iter
        self.entropy_scores = []
        self.jitter         = jitter
        self.window         = window
        self.score_vec      = score_vec
    def _first_profile_generation(self):
        """
        アルゴリズムの一番最初に呼び出され，最初のprofileを作る

        """
        simmat = np.zeros((self.nmats, self.nmats))
        # 上三角行列（対角成分除く）の要素を計算
        for i in range(self.nmats-1):
            for j in range(i+1, self.nmats):
                simmat[i, j] = ceval_edit_sim_without_jitter_virtual(self.mats[i].astype(np.float32),
                                                                     self.mats[j].astype(np.float32),
                                                                     self.score_vec.astype(np.float32))
        most_similar_pair = np.unravel_index(np.argmax(simmat), simmat.shape)
        i, j = most_similar_pair
        self.processed_ids[i] = 0
        self.processed_ids[j] = 1
        self.processed_id = 2
        self.processed[i] = True
        self.processed[j] = True
        dp, bp = ceval_edit_sim_bp_without_jitter_virtual(self.mats[i].astype(np.float32),
                                                          self.mats[j].astype(np.float32),
                                                          self.score_vec.astype(np.float32))
        self.mats[i], self.mats[j] = align_edit_sim(self.mats[i].astype(np.float32),
                                                    self.mats[j].astype(np.float32), bp.astype(np.int32))
        num_row = self.mats[i].shape[0]
        max_col = max(self.mats[i].shape[1], self.mats[j].shape[1])
        print('max_col', max_col)
        self.profile = generate_profile(num_row, max_col, False, self.mats[i], self.mats[j])
        print('profile shape', self.profile.shape)

    def _iterative_update(self):
        """
        profileに一番近いmatを用いてprofileを更新
        """
        simvec = np.zeros(self.nmats)
        unprocessed = np.where(self.processed == False)[0]
        num_row = self.mats[0].shape[0]
        if len(unprocessed) != 0:
            for i in unprocessed:
                print('shapes', self.mats[i].shape, self.profile.shape, end=' ')
                simvec[i] = ceval_edit_sim_without_jitter_virtual(self.mats[i].astype(np.float32),
                                self.profile.astype(np.float32), self.score_vec.astype(np.float32))
            next_id = np.argmax(simvec)
            self.processed[next_id] = True
            mat = self.mats[next_id]
            self.processed_ids[next_id] = self.processed_id
            self.processed_id += 1
            dp, bp = ceval_edit_sim_bp_without_jitter_virtual(mat.astype(np.float32),
                                            self.profile.astype(np.float32), self.score_vec.astype(np.float32))
            self.mats[next_id], profile = align_edit_sim(mat.astype(np.float32), self.profile.astype(np.float32), bp.astype(np.int32))
            max_col = max(self.mats[next_id].shape[1], self.profile.shape[1])
            # 無駄が多い
            tmp_mats = [mat for idx, mat in enumerate(self.mats) if self.processed[idx] == True]
            num_cols = [mat.shape[1] for mat in self.mats] + [profile.shape[1]]
            max_col = max(num_cols)
            self.profile = generate_profile(num_row, max_col, False, *tmp_mats)
        else:
            next_id = np.argmin(self.processed_ids)
            # 無駄が多い
            tmp_mats = [mat for idx, mat in enumerate(self.mats) if idx != next_id]
            num_cols = [mat.shape[1] for mat in tmp_mats] + [self.profile.shape[1]]
            max_col = max(num_cols)
            self.profile = generate_profile(num_row, max_col, True, *tmp_mats)
            mat = self.mats[next_id]
            self.processed_ids[next_id] = self.processed_id
            self.processed_id += 1
            dp, bp = ceval_edit_sim_bp_without_jitter_virtual(mat.astype(np.float32),
                                                    self.profile.astype(np.float32),
                                                    self.score_vec.astype(np.float32))
            self.mats[next_id], self.profile = align_edit_sim(mat.astype(np.float32), self.profile.astype(np.float32), bp.astype(np.int32))
            num_cols = [mat.shape[1] for mat in self.mats] + [self.profile.shape[1]]
            max_col = max(num_cols)
            self.profile = generate_profile(num_row, max_col, True, *self.mats)
            self.entropy_scores.append(cal_entropy_score(self.profile))

    def fit(self):
        self._first_profile_generation()
        for i in range(self.max_iter):
            self._iterative_update()
        self.entropy_scores = np.array(self.entropy_scores)
        # 最終的に得られたprofileと各行列のLCSを計算し，格納する．これがうまいことsequenceになってくれるはず...
        for i in range(self.nmats):
            dp, bp = ceval_edit_sim_bp_without_jitter_virtual(self.profile.astype(np.float32),
                                                                self.mats[i].astype(np.float32),
                                                    self.score_vec.astype(np.float32))
            self.lcss.append(
                cprint_LCS_without_jitter_virtual(
                    bp.astype(np.int32), self.profile.astype(np.float32), self.mats[i].astype(np.float32)))


# util関数
def generate_profile(num_row, max_col, flag, *args):
    cum_mat = np.zeros((num_row, max_col))
    shapes = [mat.shape for mat in args]
    for mat in args:
        ncol = mat.shape[1]
        cum_mat[:, :ncol] += mat
    # regularization
    row_sum = np.sum(cum_mat, axis=0)
    for col in range(cum_mat.shape[1]):
        if row_sum[col] != 0:
            cum_mat[:, col] /= row_sum[col]
    print(flag)
    if flag == True:
        cum_mat = eliminate_spaces_from_profile(cum_mat, spaces=100)
    plt.imshow(cum_mat, interpolation='nearest')
    plt.colorbar()
    plt.show()
    return cum_mat

def cal_entropy_score(cum_mat):
    # calculate entropy for each col
    E_col = np.zeros(cum_mat.shape[1])
    for col in range(cum_mat.shape[1]):
        for row in range(cum_mat.shape[0]):
            # just to avoid the calculation of log(0)
            if cum_mat[row, col] != 0:
                E_col[col] += cum_mat[row, col] * np.log(cum_mat[row, col])
    return np.sum(E_col)

def align_edit_sim(mat1, mat2, bp):
    zerovec1 = np.zeros(mat1.shape[0])
    zerovec2 = np.zeros(mat2.shape[0])
    i = bp.shape[0] - 1
    j = bp.shape[1] - 1
    mat1_list = list(mat1.T)
    mat2_list = list(mat2.T)
    alignment1 = []
    alignment2 = []
    while True:
        if i==0 or j == 0:
            break
        if bp[i, j] == 2:
            alignment1.append(mat1_list[i-1])
            alignment2.append(mat2_list[j-1])
            i -= 1
            j -= 1
        # 0 1を入れ替えた 2015/11/03 22:08
        elif bp[i, j] == 0:
            alignment1.append(mat1_list[i-1])
            alignment2.append(zerovec2)
            i -= 1
        elif bp[i, j] == 1:
            alignment1.append(zerovec1)
            alignment2.append(mat2_list[j-1])
            j -= 1
    return np.fliplr(np.array(alignment1).T), np.fliplr(np.array(alignment2).T)

def eliminate_spaces_from_profile(profile, spaces=100):
    tmat_alt = np.zeros_like(profile)
    col_alt = 0
    col = 0
    while col < profile.shape[1]:
        if np.sum(profile[:, col:(col+spaces)]) != 0:
                tmat_alt[:, col_alt] = profile[:, col]
                col_alt += 1
                col += 1
        else:
            col += spaces
    # 後処理
    col = 0
    tmat_reduced = None
    while col < profile.shape[1]:
        if np.sum(tmat_alt[:, col:]) != 0:
            col += 1
        else:
            tmat_reduced = np.zeros((profile.shape[0], col))
            tmat_reduced[:, :] = tmat_alt[:, :col]
            return tmat_reduced
    if tmat_reduced == None:
        return profile
