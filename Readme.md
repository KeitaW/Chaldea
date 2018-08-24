# Chaldea
CLI tools for sequential activity detection in neuronal data which are used for the analysis in the following paper,

Watanabe, K., Haga, T., Euston, D. R., Tatsuno, M., & Fukai, T. (2017). Unsupervised detection of cell-assembly sequences with edit similarity score. bioRxiv, 202655. http://doi.org/10.1101/202655

__This repo is deprecated.__ You can find the latest implementation in [this repo](https://github.com/KeitaW/spikesim).

Algorithms are described in [bioRxiv paper](https://www.biorxiv.org/content/early/2017/10/30/202655).

Note that the implementation of clustering algorithms (Optics and Copra) are not included in this repo. You can get the original code from the links shown below:

* Copra: http://gregory.org/research/networks/software/copra.html
* OPTICS: https://github.com/espg/OPTICS/blob/master/OPTICS.py

## Requirements
* Python 3
* Julia 0.6
* openjdk "1.8.x"

## Brief tutorial
## Store your data into __data__ directory
You need to store your data in __data__ directory. You need two different file for one experiment.
* clu: contain neuron indices. The first row is the number of neurons in it.
* res: contain spike timings unit should be (ms).
Note that if the number of data points in res file are $N$, it should be $N+1$ in clu file.
And if you have multiple recorded data from multiple shanks, you can process them simultaneously if these file names are defined as following manner,

```
all.clu.0, all.res.0, all.clu.1, all.res.1, ... all.clu.N, all.res.N
```

## Activate environment
To set PATH to the commands in __bin__, you need to execute the following line:
```bash
source activate
```
If you want to resume old environment you may run:
```bash
source deactivate
```
or just re-open your terminal.

## By the way...
You can see usage of each command with `--help` option.

## Set Root directory for your data analysis
In Chaldea, all processed data will be stored in __results__. Everytime you start analysis for a novel data set, you need initialization:
```
init_topdir ../data/tutorial
```
which will yield root directory (which should be like `../results/tutorial_170202T141337`) for successive analysis.


## Add session data into the root directory
You can add session data with `add_session` command.
```
add_session ../data/tutorial/0nneurons_10_seq_duration_100_overlap_0.0 ../results/tutorial_170202T141337
```
which will create session directory (ex. `0nneurons_10_seq_duration_100_overlap_0.0_170202T152716`) which contains the following data.

* __activity.npz__:
* __log.txt__:

## Generate binarray
In Chaldea, binned spike matrix (called binarray) will be used. Following command will create it with user-specified bin width.
```
 generate_binarray --binwidth 1 ../results/tutorial_170202T141337/0nneurons_10_seq_duration_100_overlap_0.0_170202T152716
 ```

 ## Generate simmat
 Similarity matrix (called simmat. See the article for detail) will be generated by `generate_simmat` command:

 ```
 generate_simmat ../results --a 0.05 --p 1 --window 200 --mlenseq 5
 ```
 Note: this step will take long time.

 Options are shown in below.
 * __a__: strength of exponentially growing gap penalty (optional)
 * __p__: number of parallel jobs (optional)
 * __window__: length of sliding time window(ms)
 * __mlenseq__: minimum length of sequence used for calculation time reduction.
 * __background__: run the command in the background
 * __cluster__: run the command in the cluster system

 Note: to use __cluster__ option, you need to configure a template job script in `bin/cluster` directory. A sample script is in there by default.

 ## Clustering
 The simmat has of feature space of sliding time windows that represeant relasionships of each time window. Following clustering enable us to extract time windows that are similar and repetedly occured.
 ```
 clustering_simmat ../results/tutorial_170202T141337/0nneurons_10_seq_duration_100_overlap_0.0_170202T152716/bin_size1_170202T155547/simmat_window_100a_0.5min_len_3_20170202T195657 --MinPts 5 --v 2
 ```
 * __MinPts__: Minimum criteria of cluster. Used in OPTICS.
 * __v__: Maximum number of lables that each data points retain. Used in COPRA.
 You can make clustering visualization:
 ```
 visualize_clustering ../results/tutorial_170202T141337/0nneurons_10_seq_duration_100_overlap_0.0_170202T152716/bin_size1_170202T155547/simmat_window_100a_0.5min_len_3_20170202T195657/clusters_MinPts10_v10_170203T100154
 ```
 Following command will launch small web server to see clustering resutls.
 ```
 view_clustering_results ../results/tutorial_170202T141337/0nneurons_10_seq_duration_100_overlap_0.0_170202T152716/bin_size1_170202T155547/simmat_window_100a_0.5min_len_3_20170202T195657
 ```
 Note: Currently, you may have trouble with this command when you use Windows operating system.

 ## Generate Profile
 From clustering results, common temporal structure in each cluster can be extracted by the following command:
 ```
 generate_profile ../results/tutorial_170202T141337/0nneurons_10_seq_duration_100_overlap_0.0_170202T152716/bin_size1_170202T155547/simmat_window_100a_0.5min_len_3_20170202T195657/clusters_MinPts8_v3_170203T095951 --p 3 --niter 1000
 ```

 ## Extract Sequence
 You can extract sequences by taking common temporal structure between actual spiking activity and the profiles.
 ```
 extract_sequence ../results/tutorial_170202T141337/0nneurons_10_seq_duration_100_overlap_0.0_170202T152716/bin_size1_170202T155547/simmat_window_100a_0.5min_len_3_20170202T195657/clusters_MinPts8_v3_170203T095951/profiles_numiter_1000_20170203T102604
 ```
 And you can automatically visualize the results with the following command.
 ```
 visualize_sequence ../results/tutorial_170202T141337/0nneurons_10_seq_duration_100_overlap_0.0_170202T152716/bin_size1_170202T155547/simmat_window_100a_0.5min_len_3_20170202T195657/clusters_MinPts8_v3_170203T095951/profiles_numiter_1000_20170203T102604/sequences_hosei0.0_20170203T104208/
 ```
 Also the following command launch a small web browser to check the results.

 ## Check Dtected Sequences (with some additional figures)

