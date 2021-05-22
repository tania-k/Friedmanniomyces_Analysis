#!/usr/bin/bash
#SBATCH -p short -N 1 -n 24 --mem 48gb --out logs/cds_yn00_all.log

module unload perl
module load parallel
module load subopt-kaks

INDIR=CDS_aligned_per_species
# needs biopython - so may need to module load something like funannotate
python3 scripts/gather_orthogroups_duplicates.py

OUTDIR=pairwise_ks
mkdir -p $OUTDIR
for d in $(ls $INDIR)
do
	cp lib/yn00_align_header.tsv $OUTDIR/$d.KaKs.tsv
	parallel -j 24 yn00_cds_optimal --showtable --global --noheader {} ::: $(ls $INDIR/$d/*.cds.fas)i | grep -v Stop | sed 's/[\d128-\d255]//g' >> $OUTDIR/$d.KaKs.tsv
done
