#!/bin/bash

parcores=$1

homef=/gpfs/group/asb17/default/BD2019_solar/dan
fakefilepath="fakepath/cuzidontwantpbslogfiles"

cat <<EOS | qsub -
#!/bin/bash
#PBS -l nodes=1:ppn=${parcores}
#PBS -l walltime=10:00:00
#PBS -l pmem=48gb
#PBS -A open 
#PBS -m n
#PBS -N slim$tgind 
#PBS -e "$homef"/"$fakefilepath"/logfiles/model_"$tgind".e
#PBS -o "$homef"/"$fakefilepath"/logfiles/model_"$tgind".o

cd "$homef"/solar_flare

Rscript step4_xgb_slim_tune.R "$parcores"    > \
	"$homef"/logfiles/slimtune.18-.logfile 2>&1

EOS

