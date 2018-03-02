## Prerequisite
* Java: I'm using openjdk version "1.8.0_111" but other version may also be ok.
* Python3.x:
* Julia v0.6: https://julialang.org/downloads/platform.html
This software also use bash shell scripts. So need bash also.
I'm testing this software on OSX and CentOS7.

# Dependency
This software using following two algorithm for clustering steps. They are based on existing scripts.
* OPTICS
Modified version of following two scripts are in this software (I may rewrite by myself later).
https://github.com/espg/OPTICS/blob/master/OPTICS.py
https://github.com/amyxzhang/OPTICS-Automatic-Clustering
* COPRA
The following software is used without any modification (I've just wrote a wrapper). Software http://gregory.org/research/networks/software/copra.jar
Manual
http://gregory.org/research/networks/software/copra-guide.pdf

# Installing steps
## Install Julia
You can see platform specific instructions in [here](https://julialang.org/downloads/platform.html)
Example steps for CentOS are illustrated below.
```shell
wget https://copr.fedorainfracloud.org/coprs/nalimilan/julia/repo/epel-7/nalimilan-julia-epel-7.repo
mv nalimilan-julia-epel-7.repo /etc/yum.repos.d
yum install julia
```
### Install Julia Packages
Packages required for this software are listed in `install_packages.jl` which is in top directory. These can be installed automatically with the following command.
```shell
julia install_packages.jp
```

## Install Java
 ```shell
 yum install java-1.8.0-openjdk java-1.8.0-openjdk-devel
 ```

## Install Python3
You can use [pyenv](https://github.com/pyenv/pyenv) to install various version of Python on your system.

### Install Python Packages
Packages required for this software are listed in `requirements.txt` which is in top directory. These can be installed automatically with the following command.

```shell
pip install -r requirements.txt
```

# Test library
To test behavior of the software, `test.sh` is prepared in `bin` directory.
```shell
bash test.sh
```
will execute whole analysis steps toward artificial data that is located in `data/artificial/test`.
After that, a directory that stores the result will be printed like follows:

```shell
Result is saved in  ../results/artificial_170615T203236/test_170615T203237/bin_size1_170615T203237/simmat_window_100a_0.5min_len_5_20170615T203306/clusters_MinPts50_v10_170615T203333/profiles_numiter_10_20170615T203349/sequences_sigma0_hosei0.0_20170615T203411
```

