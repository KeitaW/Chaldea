#!/bin/bash

Binwidth=1
A=0.5
WINDOW=100
P=12
MinPts=50
V=10
Niter=10

echo "Test start"
dir="./src/interim/test"; [ ! -e $dir ] && mkdir -p $dir 
echo "Preprocessing" && \
python ./src/features/preprocess_fortest.py && \
echo "Generate binarray" && \
BINARRAY_DIR=$(generate_binarray ./data/interim/test/act.npz --binwidth $Binwidth) && \
echo "Generate simmat" && \
generate_simmat $BINARRAY_DIR/binarray_data.npz --a $A --window $WINDOW --p $P && \
echo "Clustering Simmat" && \
SIMMAT_DIR=$(cat /tmp/save_dir) && \
CLUSTER_DIR=$(clustering_simmat $SIMMAT_DIR/simmat_coo.npz --MinPts $MinPts --v $V)  && \
echo "Generate profile" && \
generate_profile $CLUSTER_DIR/best-clusters-rereduced_simmat.mtx --niter $Niter  > /dev/null 2> /dev/null && \
PROFILE_DIR=$(cat /tmp/save_dir) && \
echo "Generate sequence" && \
extract_sequence $PROFILE_DIR/profiles.npz && \
SEQUENCE_DIR=$(cat /tmp/save_dir) && \
illustrate_sequences $SEQUENCE_DIR && \
echo "Test finished. See $SEQUENCE_DIR for the results"
