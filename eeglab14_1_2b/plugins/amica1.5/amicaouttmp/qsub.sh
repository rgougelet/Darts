#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe mpich 2
/home/jason/mpich2-1.5-install/bin/mpirun -np $NSLOTS -machinefile $TMP/machines /data/common/matlab/eeglab/plugins/amica1.0/amica13 /data/common/matlab/eeglab/plugins/amica1.0/amicaouttmp/input.param