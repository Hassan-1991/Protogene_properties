#Pangenome, intragenus, extragenus tagging

#Ecoli:

#Intragenus, genome:

cat Salmonella_db/Enterica_excluded_accessions.GBRS.txt Escherichia_db/Ecoli_excluded_accessions.GBRS.txt Mycobacterium_db/Tuberculosis_excluded_accessions.GBRS.txt | grep -F -f - genomes/accession_genomeID_taxonomy.tsv | sort -k2 > intragenus_genome_contig_taxa.tsv
cat Escherichia_db/Ecoli_excluded_accessions.ATB.txt Salmonella_db/Enterica_excluded_accessions.ATB.txt Mycobacterium_db/Tuberculosis_excluded_accessions.ATB.txt | grep -F -f - ATB_contigs_accessions.csv | awk -F ',' '{print $1"\t"$3"\t"$2}' >> intragenus_genome_contig_taxa.tsv

grep "Escherichia" intragenus_genome_contig_taxa.tsv | awk -F '\t' '{print $3"\t"$2}' | rev | cut -f2- -d "." | rev > Ecoli_intragenus_pangenome_contig_taxa.tsv
grep "@" /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/backup_2/all_contig_protein_taxonomy.tsv | awk -F '\t' '{print $2"\t"$1}' >> Ecoli_intragenus_pangenome_contig_taxa.tsv

sed "s/\t/,/g" Ecoli_intragenus_pangenome_contig_taxa.tsv | sort -t ',' -k2 > Ecoli_intragenus_pangenome_contig_taxa.csv


#Extragenus, genome:
sed "s/^/\t/g" odd_taxonomic_names.txt | grep -v -F -f - genomes/accession_genomeID_taxonomy.tsv | sed "s/\[//g" | sed "s/\t\'/\t/g" | sed "s/\tuncultured_/\t/g" | sed "s/\tCandidatus_/\t/g" | grep -P -v "\tbacterium" | grep -P -v "\tcandidate_" | cut -f-2 -d "_" > Ecoli_extragenus_genome_contig_taxa.tsv

#Odd taxonomic names: The "genus" designations that don't begin with upper case
alpha
arsenite-oxidising
bacterium
beta
candidate
cyanobacterium
endosymbiont
gamma
secondary
synthetic
unidentified

#ORFs:
