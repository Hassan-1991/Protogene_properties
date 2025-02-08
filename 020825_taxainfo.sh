#Pangenome, intragenus, extragenus tagging

#Ecoli:

#Intragenus, genome:

cat Salmonella_db/Enterica_excluded_accessions.GBRS.txt Escherichia_db/Ecoli_excluded_accessions.GBRS.txt Mycobacterium_db/Tuberculosis_excluded_accessions.GBRS.txt | grep -F -f - genomes/accession_genomeID_taxonomy.tsv | sort -k2 > intragenus_genome_contig_taxa.tsv
cat Escherichia_db/Ecoli_excluded_accessions.ATB.txt Salmonella_db/Enterica_excluded_accessions.ATB.txt Mycobacterium_db/Tuberculosis_excluded_accessions.ATB.txt | grep -F -f - ATB_contigs_accessions.csv | awk -F ',' '{print $1"\t"$3"\t"$2}' >> intragenus_genome_contig_taxa.tsv

grep "Escherichia" intragenus_genome_contig_taxa.tsv | awk -F '\t' '{print $3"\t"$2}' | rev | cut -f2- -d "." | rev > Ecoli_intragenus_pangenome_contig_taxa.tsv

grep "sequence-region" /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/*gff | grep -v "all_500_gffs" | cut -f-2 -d " " | rev | cut -f1 -d '/' | rev | sed "s/.gff:##sequence-region//g" > Ecoli_pangenome_contig_taxa.interim
ls /stor/work/Ochman/hassan/Ecoli_pangenome/500_gffs/*gff | grep -v "all_" | rev | cut -f1 -d '/' | rev | sed "s/.gff//g" | awk '{print $0,$0}' >> Ecoli_pangenome_contig_taxa.interim
sort -k1 Ecoli_pangenome_contig_taxa.interim -o Ecoli_pangenome_contig_taxa.interim
tail -n+2 /stor/work/Ochman/hassan/Ecoli_pangenome/500_ipp_lineagedesignations.tsv | cut -f1,5 | sort -k1 | join -1 1 -2 1 - Ecoli_pangenome_contig_taxa.interim | awk '{print $1"@"$2,$3}' >> Ecoli_intragenus_pangenome_contig_taxa.tsv
sed -i "s/Escherichia /Escherichia_/g" Ecoli_intragenus_pangenome_contig_taxa.tsv
sed -i "s/ /\t/g" Ecoli_intragenus_pangenome_contig_taxa.tsv
sort -k2 Ecoli_intragenus_pangenome_contig_taxa.tsv -o Ecoli_intragenus_pangenome_contig_taxa.tsv

#The above works for protein designations. For genome designations (i.e., found in intervalinfo file), complete contig names are required.
awk -F '\t' '{print $3"\t"$2}' Ecoli_extragenus_genome_contig_taxa.tsv > Ecoli_contig_taxonomy.intervalinfo.tsv
awk -F '\t' '{print $3"\t"$2}' Ecoli_intragenus_genome_contig_taxa.tsv >> Ecoli_contig_taxonomy.intervalinfo.tsv
tail -n+2 /stor/work/Ochman/hassan/Ecoli_pangenome/500_ipp_lineagedesignations.tsv | cut -f1,5 | sort -k1 | join -1 1 -2 1 - Ecoli_pangenome_contig_taxa.interim | awk '{print $1"@"$2,$3}' >> Ecoli_contig_taxonomy.intervalinfo.tsv
sed -i "s/Escherichia /Escherichia_/g" Ecoli_contig_taxonomy.intervalinfo.tsv
sed -i "s/ /\t/g" Ecoli_contig_taxonomy.intervalinfo.tsv
sort -k2 Ecoli_contig_taxonomy.intervalinfo.tsv -o Ecoli_contig_taxonomy.intervalinfo.tsv


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
