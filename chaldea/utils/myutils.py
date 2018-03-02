import numpy as np
import os
from scipy.sparse import coo_matrix
import re
import json
import xml.etree.ElementTree as etree
import pandas as pd

def ret_dotdot(path, up=1):
    original_path = os.getcwd()
    os.chdir(path)
    for i in range(up):
        os.chdir("..")
    ret_dir = os.getcwd()
    os.chdir(original_path)
    return ret_dir

def load_copra_data(rf, drop=5):
    data = []
    for line in open(rf,'r').readlines():
        tmp = [int(l) for l in line[:-2].split(' ')]
        if len(tmp) >= drop:
            data.append([int(l)-1 for l in line[:-2].split(' ')])
    return data

def ret_reference_points_list(correspondence_table_time, data, binsize):
    reference_points_list = []
    for d in data:
        reference_points_list.append(correspondence_table_time[d] / (10**int(3-np.log10(binsize))))
    reference_points_list = np.array(reference_points_list)
    return reference_points_list

def ret_coo_matrix_from_data(data):
    return coo_matrix((data["data"], (data["row"], data["col"])))

def ret_data_from_clust(clusters_dir):
    # read data from clusters_dir
    rf = os.path.join(clusters_dir, "best-clusters-rereduced_simmat.mtx")
    data = load_copra_data(rf, drop=5)
    correspondence_table_time = np.load(os.path.join(clusters_dir, "correspondence_table_time.npz"))["arr_0"] - 1
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
    return (data, correspondence_table_time, reference_points_list, simmat_data, params, binarray_data, binsize)

def ret_data(profiles_dir):
    # read data from clusters_dir
    profiles = np.load(os.path.join(profiles_dir, "profiles.npz"))
    clusters_dir = ret_dotdot(profiles_dir)
    rf = os.path.join(clusters_dir, "best-clusters-rereduced_simmat.mtx")
    data = load_copra_data(rf, drop=5)
    correspondence_table_time = np.load(os.path.join(clusters_dir, "correspondence_table_time.npz"))["arr_0"] - 1
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
    act = np.load(os.path.join(session_dir, "act.npz"))
    return (profiles, data, correspondence_table_time, reference_points_list, simmat_data, params, binarray_data)

def load_sequence_dicts(sequence_rows_dict, sequence_cols_dict):
    '''
    sequence_rows, sequence_colsはdictionary formatで保存されているがそれをそのまま使うのはちょっと不便．
    なのでリストの形に変換してやる．
    '''
    all_keys = sequence_cols_dict.keys()
    num_seq_kinds = max([int(key.split('.', 1)[0]) for key in all_keys])
    sequence_rows_list = []
    sequence_cols_list = []
    for seq_idx in range(1, num_seq_kinds+1):
        matches = [re.findall(r'^'+str(seq_idx)+r'\.[0-9]+', key) for key in all_keys]
        seq_indices = [key for key in matches if len(key) != 0]
        sequence_cols = [sequence_cols_dict[key[0]] for key in seq_indices]
        sequence_rows = [sequence_rows_dict[key[0]] for key in seq_indices]
        sequence_rows_list.append(sequence_rows)
        sequence_cols_list.append(sequence_cols)
    return sequence_rows_list, sequence_cols_list

# ちょっとスマートなやりかたと違うけれども...profileまでではなくsequenceまでをいっきに読み込む関数
def ret_sequences(sequence_dir, actfile="act.npz", drop=5):
    # read data from sequence_dir
    sequence_rows_dict = np.load(os.path.join(sequence_dir, "sequence_rows.npz"))
    sequence_cols_dict = np.load(os.path.join(sequence_dir, "sequence_cols.npz"))
    sequence_rows, sequence_cols = load_sequence_dicts(sequence_rows_dict, sequence_cols_dict)
    profiles_dir = ret_dotdot(sequence_dir)
    # read data from profiles_dir
    profiles = np.load(os.path.join(profiles_dir, "profiles.npz"))
    # read data from clusters_dir
    clusters_dir = ret_dotdot(profiles_dir)
    rf = os.path.join(clusters_dir, "best-clusters-rereduced_simmat.mtx")
    data = load_copra_data(rf, drop=drop)
    correspondence_table_time = np.load(os.path.join(clusters_dir, "correspondence_table_time.npz"))["arr_0"] - 1
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
    act = np.load(os.path.join(session_dir, actfile))
    return (sequence_rows, sequence_cols, profiles, data, correspondence_table_time, reference_points_list, simmat_data, params, binarray_data)

# ココらへんには，CRCNSのLiner maze taskで取れたデータに対してのみ有効な関数がちらほら見受けられる
# そういったものは別のモジュールに切り分けたほうがいいかも知らん

def decompose_mixed_eeg_data(eeg_mixed, num_channels):
    '''
    In *.eeg file... 説明書きを追加すること
    また，今後も使いそうだから，util関数に追加すること
    '''
    nrow = eeg_mixed.shape[0] // num_channels
    ncol = num_channels
    eeg_separated = np.zeros((nrow, ncol), dtype = eeg_mixed.dtype)
    for channel in range(num_channels):
        eeg_separated[:, channel] = eeg_mixed[channel::num_channels]
    return eeg_separated

def read_eeg(datadir, num_channels, sample_rate):
    eeg_file = os.path.join(datadir, datadir.split('/')[-1] + ".eeg")
    eeg_mixed = np.fromfile(eeg_file, dtype='<i2') / 1000
    eeg_separated = decompose_mixed_eeg_data(eeg_mixed, num_channels)
    eeg_time = np.arange(eeg_separated.shape[0]) / sample_rate
    return eeg_time, eeg_separated

def acquire_nChannels(datadir):
    xml_file = os.path.join(datadir, datadir.split('/')[-1] + ".xml")
    tree = etree.parse(xml_file)
    root = tree.getroot()
    acqisitionSystem = root.findall("acquisitionSystem")[0]
    nChannels = acqisitionSystem.findall("nChannels")[0]
    for it in nChannels.itertext():
        num_channels = it
    return int(num_channels)

def acquire_eeg_sample_rate(datadir):
    xml_file = os.path.join(datadir, datadir.split('/')[-1] + ".xml")
    tree = etree.parse(xml_file)
    root = tree.getroot()
    fieldPotentials = root.findall("fieldPotentials")[0]
    lfpSamplingRate = fieldPotentials.findall("lfpSamplingRate")[0]
    for lsr in lfpSamplingRate.itertext():
        sample_rate = lsr
    return int(sample_rate)

def acquire_shank_from_sequence_dir(sequences_dir):
    binarray_dir = ret_dotdot(sequences_dir, 4)
    binarray_dir = binarray_dir.split("/")[-1]
    shank = re.findall(r'_e[0-9]+_', binarray_dir)[0]
    shank = shank[1:-1]
    return shank

def acquire_region_from_BasicInfo(sequences_dir, shank):
    session_dir = ret_dotdot(sequences_dir, 5)
    BasicInfo = pd.read_csv(os.path.join(session_dir, "BasicInfo.tsv"), sep="\t")
    return(BasicInfo[shank][0])

def acquire_channels_corresponding_to_shank(data_dir, shank):
    with open(os.path.join(data_dir, "channel_info.json"), 'r') as rf:
        channel_info = json.load(rf)
    return [int(ci) for ci in channel_info[shank[1:]]]

def get_data_dir_from_result_directory(result_dir):
    tmp = re.findall(r'ec[0-9]+.[0-9]+', result_dir)
    data_dir = "../data/crcns/" + tmp[0] + "/" + tmp[1]
    return data_dir

def read_whl(datadir):
    whlfile = re.split(r"\/", datadir)[-1] + ".whl"
    #act = np.load(os.path.join(activity_dir, "activity.npz"))
    whl = np.loadtxt(os.path.join(datadir, whlfile))
    whl = pd.DataFrame(whl)
    whl.columns = ["x1", "y1", "x2", "y2"]
    whl[whl == -1] = np.NaN
    x1 = whl["x1"]
    x2 = whl["x2"]
    times = np.linspace(0, len(x1)/39.06, len(x1))
    return times, x1, x2

def ret_datadir(clusters_dir):
    extract_paths = (lambda dir_: (lambda paths: (paths[2], paths[3]))(re.split(r"\/", dir_)))
    remove_date = (lambda dir_: (lambda paths: "_".join(paths))(re.split(r"_", dir_)[:-1]))
    topdir, sessiondir = [remove_date(dirname) for dirname in extract_paths(clusters_dir)]
    datadir = os.path.join("..", "data", "crcns", topdir, sessiondir)
    return(datadir)

def ret_puseudo_cludet(clusters_dir):
    binarray_dir = ret_dotdot(clusters_dir, 2)
    binarray_data = np.load(os.path.join(binarray_dir, "binarray_data.npz"))
    binsize = int(binarray_data["binsize(ms)"])
    rf = os.path.join(clusters_dir, "best-clusters-rereduced_simmat.mtx")
    data = load_copra_data(rf, drop=80)
    correspondence_table_time = np.load(os.path.join(clusters_dir, "correspondence_table_time.npz"))["arr_0"] - 1
    clu = list(); det = list()
    for idx, d in enumerate(data):
        part_det = list(correspondence_table_time[np.array(d)] / (10**int(3-np.log10(binsize))))
        part_clu = list(np.ones(len(part_det)) * idx)
        clu += part_clu; det += part_det
    clu = np.array(clu); det = np.array(det)
    sorted_indices = np.argsort(det)
    return clu[sorted_indices], det[sorted_indices]
