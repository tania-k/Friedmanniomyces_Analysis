#!/usr/bin/bash
#SBATCH -N 1 -n 8 --mem 2gb -p short --out logs/makepairs.%a.v3.log -J makepairs

module unload perl
module load perl/5.22.0

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
	CPUS=1
fi

PEP=Aligned_PEP4
CDS=results4
Al_CDS=Aligned_CDS4
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
    FILE=$PEP/$NAME.fa
    echo $FILE
    if [ ! -f $Al_CDS/$NAME.cds.aln ]; then
        ./scripts/bp_mrtrans.pl -i $FILE -if fasta -s $CDS/$NAME.CDS.fa.new -of fasta -o $Al_CDS/$NAME.cds.aln
    fi
done
