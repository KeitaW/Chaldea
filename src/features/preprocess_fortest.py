import subprocess
import os
import sys
import json
import numpy as np

chaldeadir  = os.environ.get("CHALDEADIR")
dsessiondir = os.path.join(chaldeadir, "data/raw/test")
clufile     = os.path.join(chaldeadir, "data/raw/test/all.clu.1")
resfile     = os.path.join(chaldeadir, "data/raw/test/all.res.1")
rsessiondir = os.path.join(chaldeadir, "data/interim/test")
params_file = os.path.join(chaldeadir, "data/raw/test/params.json")
if not os.path.exists(rsessiondir):
    os.makedirs(rsessiondir)
with open(params_file, "r") as rf:
    params = json.load(rf)

duration_s = params["duration(ms)"] / 1000   
nneurons   = params["nneurons"]
clu = np.loadtxt(clufile)
res = np.loadtxt(resfile) / 1000

# * neuronids: neuronid
# * spikes: time of spike (s)
# * nneuorns: total number of neurons
# * duration_s: duration (s)
savedata = {"duration_s": duration_s, "nneurons": nneurons, "neurons": clu[1:], "spikes": res}
np.savez(os.path.join(rsessiondir, "act.npz"), **savedata)
