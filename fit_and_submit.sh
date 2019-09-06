#!/bin/bash

parcores=$1
tgind=$2

homef=/gpfs/group/asb17/default/BD2019_solar/dan
fakefilepath="fakepath/cuzidontwantpbslogfiles"

cat <<EOS | qsub -
#!/bin/bash
#PBS -l nodes=1:ppn=${parcores}
#PBS -l walltime=2:00:00
#PBS -l pmem=48gb
#PBS -A open 
#PBS -m n
#PBS -N modno$tgind 
#PBS -e "$homef"/"$fakefilepath"/logfiles/model_"$tgind".e
#PBS -o "$homef"/"$fakefilepath"/logfiles/model_"$tgind".o

cd "$homef"/solar_flare

Rscript fit_and_submit.R "$parcores" "$tgind"   > \
	"$homef"/logfiles/fitsubmit."$tgind".logfile 2>&1

EOS

