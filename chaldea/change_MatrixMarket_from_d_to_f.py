#-*-coding:utf-8-*-
from __future__ import division
import numpy as np
import sys
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
import os
# juliaのmmread関数は残念なことにFloat64でないと読めないっぽい．Int型で作ってしまった行列をFloat64に変換してやるスクリプト
filename = sys.argv[1]
binarray_coo = mmread(filename)
rows = binarray_coo.row
cols = binarray_coo.col
datas = binarray_coo.data
source_dir = os.path.dirname(filename)
mmwrite(filename, binarray_coo, field='real')
np.save(os.path.join(source_dir+'rows.npz'), rows)
np.save(os.path.join(source_dir+'cols.npz'), cols)
np.save(os.path.join(source_dir+'datas.npz'), datas)
