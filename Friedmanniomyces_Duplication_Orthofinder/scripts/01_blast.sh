#!/bin/bash
#SBATCH  --time 1-0:00:00 --ntasks 8 --nodes 1 --mem 24G --out logs/fasta.v4.%a.log

module load ncbi-blast/2.2.30+

CPU=1
if [ $SLURM_CPUS_ON_NODE ] ; then
    CPU=$SLURM_CPUS_ON_NODE
fi


INFILE=new.txt
N=${SLURM_ARRAY_TASK_ID}
INDIR=seqIDs4
OUT=results4

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
    FILE=$INDIR/$NAME.fa.seqID.cleaned.seqID
    echo $FILE
    if [ ! -f $OUT/$NAME.CDS.fa ]; then
        blastdbcmd -entry_batch "$FILE" -db Dothid >> $OUT/$NAME.CDS.fa 
    fi
done
