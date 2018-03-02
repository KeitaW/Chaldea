#!/bin/bash
# this is a example script that makes generate_simmat command runs in a cluster node.

# https://stackoverflow.com/questions/3572030/bash-script-absolute-path-with-osx
# as usual SCRIPT_DIR=$(cd $(dirname $0); pwd) not working correctly in OSX (it contains PROMPT for some reason...)
realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}


BINARRAY=$1
WINDOW=$2
A=$3
SLIDE=$4
P=$5
MLENSEQ=$6
EXFLAG=$7
HOST=$8
SCRIPT="jobscript.sh"

cat <<EOT > $SCRIPT
#!/bin/bash
#$ -cwd
#$ -q $HOST
#$ -M 092975@gmail.com
#$ -m abe
#$ -pe OpenMP $P

source ~/.bashrc
source $CHALDEA/bin/activate
cd $CHALDEA/bin
echo "julia -p $P $CHALDEA/bin/pmap_cal_simmat.jl $BINARRAY $A $WINDOW $SLIDE $MLENSEQ $EXFLAG"
julia -p $P $CHALDEA/bin/pmap_cal_simmat.jl $BINARRAY $A $WINDOW $SLIDE $MLENSEQ $EXFLAG

EOT

echo "The job that will be submitted..."
cat $SCRIPT
qsub $SCRIPT
rm $SCRIPT
