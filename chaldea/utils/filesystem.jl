using NPZ
include("load_copra_data.jl")

function ret_dotdot(path)
    original_path = pwd()
    cd(path)
    cd("..")
    ret_dir = pwd()
    cd(original_path)
    return ret_dir
end
 
function ret_dotdot(path, up=1)
    original_path = pwd()
    cd(path)
    for i in 1:up
        cd("..")
    end
    ret_dir = pwd()
    cd(original_path)
    return ret_dir
end

function ret_reference_points_list(correspondence_table_time, data, binsize)
    reference_points_list = []
    for d in data
        push!(reference_points_list, correspondence_table_time[d] / (10^Int64(3-log10(binsize))))
    end
    return reference_points_list
end

function ret_coo_matrix_from_data(data)
    return sparse(data["row"],  data["col"], data["data"])
end

function ret_data(profiles_dir)
    # read data from clusters_dir
    profiles = npzread(joinpath(profiles_dir, "profiles.npz"))
    clusters_dir = ret_dotdot(profiles_dir)
    rf = joinpath(clusters_dir, "best-clusters-rereduced_simmat.mtx")
    data = load_copra_data(rf, drop=5)
    correspondence_table_time = npzread(joinpath(clusters_dir, "correspondence_table_time.npz"))["arr_0"] - 1
    # read data from simmat_dir    
    simmat_dir = ret_dotdot(clusters_dir)
    simmat_data = npzread(joinpath(simmat_dir, "simmat_coo.npz"))
    simmat_data["row"] += 1
    simmat_data["col"] += 1
    params = npzread(joinpath(simmat_dir, "parameters.npz"))
    # read data from binarray_dir        
    binarray_dir = ret_dotdot(simmat_dir)
    binarray_data = npzread(joinpath(binarray_dir, "binarray_data.npz"))
	binarray_data["row"] += 1
	binarray_data["col"] += 1
    binsize = Int64(binarray_data["binsize(ms)"])
    reference_points_list = ret_reference_points_list(correspondence_table_time, data, binsize)
    # read data from session_dir
    session_dir = ret_dotdot(binarray_dir)
    #act = npzread(joinpath(session_dir, "activity.npz"))
    return (profiles, data, correspondence_table_time, binsize, reference_points_list, simmat_data, params, binarray_data)
end

# ちょっとスマートなやりかたと違うけれども...profileまでではなくsequenceまでをいっきに読み込む関数
function ret_sequences(sequence_dir)
    # read data from sequence_dir
    sequence_rows = npzread(joinpath(sequence_dir, "sequence_rows.npz"))
    sequence_cols = npzread(joinpath(sequence_dir, "sequence_cols.npz"))
    profiles_dir = ret_dotdot(sequence_dir)
    # read data from profiles_dir
    profiles = npzread(joinpath(profiles_dir, "profiles.npz"))
    # read data from clusters_dir
    clusters_dir = ret_dotdot(profiles_dir)
    rf = joinpath(clusters_dir, "best-clusters-rereduced_simmat.mtx")
    data = load_copra_data(rf, drop=5)
    correspondence_table_time = npzread(joinpath(clusters_dir, "correspondence_table_time.npz"))["arr_0"] - 1
    # read data from simmat_dir    
    simmat_dir = ret_dotdot(clusters_dir)
    simmat_data = npzread(joinpath(simmat_dir, "simmat_coo.npz"))
    params = npzread(joinpath(simmat_dir, "parameters.npz"))
    # read data from binarray_dir        
    binarray_dir = ret_dotdot(simmat_dir)
    binarray_data = npzread(joinpath(binarray_dir, "binarray_data.npz"))
	binarray_data["row"] += 1
	binarray_data["col"] += 1
    binsize = int(binarray_data["binsize(ms)"])
    reference_points_list = ret_reference_points_list(correspondence_table_time, data, binsize)
    # read data from session_dir
    session_dir = ret_dotdot(binarray_dir)
    act = npzread(joinpath(session_dir, "activity.npz"))
    return (sequence_rows, sequence_cols, profiles, data, correspondence_table_time, reference_points_list, simmat_data, params, binarray_data)
end

