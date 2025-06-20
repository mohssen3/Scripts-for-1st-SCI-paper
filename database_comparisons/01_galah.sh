
ml Galah/0.3.1

galah cluster --genome-fasta-directory ./$input -x fa --checkm-tab-table ./${input}_checkm.txt\
 --ani 95 --precluster-ani 90 -t 96\
 --output-cluster-definition ./${input}_representative_genomes.tsv\
 --output-representative-fasta-directory-copy ./${input}_representative_genomes

## change --ani 95 to --ani 99 if strain-level