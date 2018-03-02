import numpy as np
import scipy as sp
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_style({"grid.linestyle": ""})
sns.set_context("paper", font_scale=2)
import re
import datetime
import scipy.fftpack
from scipy import signal
import json
from pylab import get_cmap
from scipy import stats
from itertools import chain
import pandas as pd
from itertools import cycle
from scipy import io
import matplotlib.gridspec as gridspec
from sklearn.manifold import Isomap
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os, sys
from sklearn.decomposition import PCA, FastICA
sys.path.append("../lib/utils/")
import myutils
sys.path.append("../lib")
from f_value import f_value
import csv
from sklearn.metrics import v_measure_score, adjusted_mutual_info_score

# Data dir
#clusters_dir = "../results/artificial_170307T204910/test_assembly_170320T162851/bin_size10_170320T162920/simmat_window_100a_0.8_20170320T165040/clusters_MinPts2_v2_170320T165250"
def ret_datadir(cluster_dir):
    act_dir = myutils.ret_dotdot(cluster_dir, 3)
    data_dir = ""
    for path in re.split("_", os.path.basename(act_dir))[:-1]:
        data_dir += path
        data_dir += "_"
    data_dir = os.path.join("../data/artificial", data_dir[:-1])
    return(data_dir)

def ret_dotdot(path, up=1):
    original_path = os.getcwd()
    os.chdir(path)
    for i in range(up):
        os.chdir("..")
    ret_dir = os.getcwd()
    os.chdir(original_path)
    return ret_dir


def load_copra_data(rf, drop=3):
    data = []
    for line in open(rf,'r').readlines():
        tmp = [int(l) for l in line[:-2].split(' ')]
        if len(tmp) >= drop:
            data.append([int(l)-1 for l in line[:-2].split(' ')])
    return data

def ret_reference_points_list(correspondence_table_time, data, binsize):
    reference_points_list = []
    for d in data:
        reference_points_list.append(correspondence_table_time[np.array(d)] / (10**int(3-np.log10(binsize))))
    reference_points_list = np.array(reference_points_list)
    return reference_points_list

def ret_coo_matrix_from_data(data):
    return coo_matrix((data["data"], (data["row"], data["col"])))

def ret_data(clusters_dir):
    # read data from clusters_dir
    rf = os.path.join(clusters_dir, "best-clusters-rereduced_simmat.mtx")
    data = load_copra_data(rf, drop=3)
    correspondence_table_time = np.load(os.path.join(clusters_dir, "correspondence_table_time.npz"))["arr_0"] 
    # read data from simmat_dir    
    simmat_dir = ret_dotdot(clusters_dir)
    simmat_data = np.load(os.path.join(simmat_dir, "simmat_coo.npz"))
    params = np.load(os.path.join(simmat_dir, "parameters.npz"))
    # read data from binarray_dir        
    binarray_dir = ret_dotdot(simmat_dir)
    binarray_data = np.load(os.path.join(binarray_dir, "binarray_data.npz"))
    binsize = int(binarray_data["binsize(ms)"])
    reference_points_list = ret_reference_points_list(correspondence_table_time, data, binsize)
    # read data from session_dir
    session_dir = ret_dotdot(binarray_dir)
    act = np.load(os.path.join(session_dir, "activity.npz"))
    return (data, correspondence_table_time, reference_points_list, simmat_data, params, binarray_data, binsize)

def take_correspondence(labels_true, labels_pred):
    """
    -1:(noise)も扱えるように改良したバージョン
    """
    clust_labels = np.unique(labels_true)
    corresponding_labels_pred = np.ones_like(labels_pred) * (-1)
    correspondence_table = dict()

    for label in clust_labels[clust_labels != -1]:
        predicted_labels, labels_counts = np.unique(labels_pred[labels_true == label], return_counts=True)
        correspondence_table[label] = predicted_labels[np.argmin(np.abs(labels_counts-24))]
        #correspondence_table[label] = predicted_labels[np.argmax(labels_counts)]
        corresponding_labels_pred[labels_pred == correspondence_table[label]] = label
    for label in np.unique(labels_pred):
        if label in correspondence_table.values():
            continue
        else:
            corresponding_labels_pred[labels_pred == label] = -1            
    return(corresponding_labels_pred)
import json
def load_org_data(data_dir):
    with open(os.path.join(data_dir, "params.json"), "r") as rf:
        params = json.load(rf)
    return(np.loadtxt(dtype=np.int, fname=os.path.join(data_dir, "all.clu.1")),
           np.loadtxt(dtype=np.float, fname=os.path.join(data_dir, "all.res.1")),
           np.loadtxt(dtype=np.int, fname=os.path.join(data_dir, "all.idx.1")),
           np.loadtxt(dtype=np.int, fname=os.path.join(data_dir, "onset_of_sequences.txt")),
           np.loadtxt(dtype=np.int, fname=os.path.join(data_dir, "onset_labels.txt")),
           params)

def ret_labels_from_clusterdir(clusters_dir):
    data_dir = ret_datadir(clusters_dir)
    clu, res, spike_indices, onset_of_sequences, onset_labels, params = load_org_data(data_dir)
    (data, correspondence_table_time, reference_points_list, simmat_data, params_, binarray_data, binsize) = ret_data(clusters_dir)
    binarray_coo = myutils.ret_coo_matrix_from_data(binarray_data)
    binarray_csc = binarray_coo.tocsc()
    binarray_csr = binarray_csc.tocsr()
    activity_dir = myutils.ret_dotdot(clusters_dir, 3)
    window_num = int(params["sequence_duration(ms)"]) // binsize
    window_ms = int(params["sequence_duration(ms)"])
    window_s = int(params["sequence_duration(ms)"]) / 1000
    simmat_dir = myutils.ret_dotdot(clusters_dir, 1)
    times = np.load(os.path.join(simmat_dir, "times.npz"))
    labels = np.ones(len(times))
    labels_true = []
    for t in times:
        if t in onset_of_sequences:
            labels_true.append(onset_labels[onset_of_sequences==t][0])
        else:
            labels_true.append(-1)
    labels_true = np.array(labels_true)

    rf = os.path.join(clusters_dir, "best-clusters-rereduced_simmat.mtx")
    clust_labels_list = load_copra_data(rf)
    labels_pred = np.ones_like(labels_true)*(-1)
    #labels_pred = np.ones(binarray_csc.shape[1] // window_num + 1) * (-1)
    for clust_idx, clust_labels in enumerate(clust_labels_list):
        for label in clust_labels:
            idx = correspondence_table_time[label]//window_num
            if int(idx) <= len(labels_pred) -1:
                labels_pred[int(idx)] = clust_idx
    return(labels_true, take_correspondence(labels_true, labels_pred))

