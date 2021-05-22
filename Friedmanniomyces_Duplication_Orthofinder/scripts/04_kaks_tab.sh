#!/usr/bin/bash
#SBATCH -N 1 -n 8 --mem 2gb -p short --out logs/kaks.%a.log -J makepairs

module unload perl
module load subopt-kaks

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
	CPUS=1
fi

PEP=Aligned_PEP4
CDS=results4
Al_CDS=Aligned_CDS4
OUTKAKS=kaks4
YN00=$(which yn00_cds_prealigned)
if [ ! -f $YN00 ]; then
    echo "need to have installed yn00_cds_prealigned - see https://github.com/hyphaltip/subopt-kaks"
    exit
fi    
INFILE=new.txt
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi

IFS=,
sed -n ${N}p $INFILE | while read NAME
do
    cp lib/yn00_header.tsv $OUTKAKS/$NAME.yn00.tab
        $YN00 $Al_CDS/$NAME.cds.aln | grep -v "SEQ1" >> $OUTKAKS/$NAME.yn00.tab
done
