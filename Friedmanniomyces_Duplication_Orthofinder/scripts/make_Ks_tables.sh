
grep -h -P '^BTJ68.*\s+BTJ68' kaks4/*.tab | cut -f4 > pairwise_ks4/Hwer.tsv
grep -h -P '^B0A54.*\s+B0A54' kaks4/*.tab | cut -f4 > pairwise_ks4/Fendolithicus.tsv
grep -h -P '^B0A55.*\s+B0A55' kaks4/*.tab | cut -f4 > pairwise_ks4/Fsimplex.tsv
