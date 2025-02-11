#Starting with the query files generated in protogene_curation

Ecoli_protein_queryfile.faa
Salmonella_protein_queryfile.faa
Mycobacterium_protein_queryfile.faa

###ECOLI###

#Outside genus, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Ecoli_vs_GBRS_annotated.tsv -k 0 -b8 -c1
#Outside genus, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Ecoli_vs_GBRS_ORFs.tsv -k 0 -b8 -c1
#Outside species, annotated 
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Ecoli_vs_noncoliEscherichia_annotated.tsv -k 0 -b8 -c1
#Outside species, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Ecoli_vs_noncoliEscherichia_ORFs.tsv -k 0 -b8 -c1

###Salmonella###

#Outside genus, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Salmonella_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_excluded.proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Salmonella_vs_GBRS_annotated.tsv -k 0 -b8 -c1
#Outside genus, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Salmonella_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_excluded.genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Salmonella_vs_GBRS_ORFs.tsv -k 0 -b8 -c1
#Outside species, annotated 
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Salmonella_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Enterica_excluded.proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Salmonella_vs_nonentericaSalmonella_annotated.tsv -k 0 -b8 -c1
#Outside species, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Salmonella_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Enterica_excluded.genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Salmonella_vs_nonentericaSalmonella_ORFs.tsv -k 0 -b8 -c1

###Mycobacterium###

#Outside genus, annotated
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Mycobacterium_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Mycobacterium_excluded.proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Mycobacterium_vs_GBRS_annotated.tsv -k 0 -b8 -c1
#Outside genus, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Mycobacterium_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Mycobacterium_excluded.genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Mycobacterium_vs_GBRS_ORFs.tsv -k 0 -b8 -c1
#Outside species, annotated 
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Mycobacterium_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Tuberculosis_excluded.proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Mycobacterium_vs_nontuberculosisMycobacterium_annotated.tsv -k 0 -b8 -c1

#Outside species, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Mycobacterium_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Tuberculosis_excluded.genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Mycobacterium_vs_nontuberculosisMycobacterium_ORFs.tsv -k 0 -b8 -c1

#pangenomes
#Ecoli:
/stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_proteins
/stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes.getorf.ATG_TTG_GTG.prot
#Salmonella:
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/all_proteins.faa --db Salmonella_pangenome_proteins
cat /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/genomes/*_genomic.fna > /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/Salmonella_pangenome.genomes.faa
makeblastdb -in /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/Salmonella_pangenome.genomes.faa -dbtype nucl -out Salmonella_pangenome.genomes
getorf -sequence /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/Salmonella_pangenome.genomes.faa -outseq /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_pangenome.getorf -table 1 -minsize 30 -find 3
seqkit fx2tab Salmonella_pangenome.getorf | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Salmonella_pangenome.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Salmonella_pangenome.genomes.getorf.ATG_TTG_GTG Salmonella_pangenome.genomes.getorf.ATG_TTG_GTG.prot.faa
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Salmonella_pangenome.genomes.getorf.ATG_TTG_GTG.prot.faa --db Salmonella_pangenome.genomes.getorf.ATG_TTG_GTG.prot
#Mycobacterium:
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_queries/Mycobacterium_331/all_331_proteins.faa --db Mycobacterium_pangenome_proteins
cat /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_queries/Mycobacterium_331/genomes/*_genomic.fna > /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_queries/Mycobacterium_331/Mycobacterium_pangenome.genomes.faa
makeblastdb -in /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_queries/Mycobacterium_331/Mycobacterium_pangenome.genomes.faa -dbtype nucl -out Mycobacterium_pangenome.genomes
getorf -sequence /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_queries/Mycobacterium_331/Mycobacterium_pangenome.genomes.faa -outseq /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_queries/Mycobacterium_331/Mycobacterium_pangenome.getorf -table 1 -minsize 30 -find 3
seqkit fx2tab /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_queries/Mycobacterium_331/Mycobacterium_pangenome.getorf | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Mycobacterium_pangenome.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Mycobacterium_pangenome.genomes.getorf.ATG_TTG_GTG Mycobacterium_pangenome.genomes.getorf.ATG_TTG_GTG.prot.faa
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in Mycobacterium_pangenome.genomes.getorf.ATG_TTG_GTG.prot.faa --db Mycobacterium_pangenome.genomes.getorf.ATG_TTG_GTG.prot

#Pangenome searches
#Ecoli
#Across pangenome, annotated
#/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_annotated.tsv -k 0 -b8 -c1
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Ecoli_vs_pangenome_annotated.tsv -k 0 -b8 -c1
#Across pangenome, ORFs
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_ORFs.tsv -k 0 -b8 -c1
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Ecoli_vs_pangenome_ORFs.tsv -k 0 -b8 -c1
#Salmonella
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Salmonella_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_pangenome_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Salmonella_vs_pangenome_annotated.tsv -k 0 -b8 -c1
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Salmonella_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_pangenome.genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Salmonella_vs_pangenome_ORFs.tsv -k 0 -b8 -c1
#Mycobacterium
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Mycobacterium_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Mycobacterium_pangenome_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Mycobacterium_vs_pangenome_annotated.tsv -k 0 -b8 -c1
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Mycobacterium_protein_queryfile.faa -d /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Mycobacterium_pangenome.genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out Mycobacterium_vs_pangenome_ORFs.tsv -k 0 -b8 -c1

for i in Mycobacterium_*
do
sed -i "s/BALROG/balrog/g" $i
sed -i "s/GMS2/gms2/g" $i
sed -i "s/PRODIGAL/prodigal/g" $i
sed -i "s/SMORFER/smorfer/g" $i
done

for i in Ecoli Salmonella Mycobacterium
do
cat "$i"_vs_GBRS_annotated.tsv "$i"_vs_GBRS_ORFs.tsv | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1 | sort -u > "$i"_step1_genusspecific_nonORFan.txt
grep -vf "$i"_step1_genusspecific_nonORFan.txt "$i"_protein_queryfile.faa | grep "^>" | tr -d ">" > "$i"_step1_genusspecific_ORFan.txt
faSomeRecords "$i"_protein_queryfile.faa "$i"_step1_genusspecific_ORFan.txt "$i"_step1_genusspecific_ORFan.faa
done

for i in Ecoli Salmonella Mycobacterium
do
cat "$i"_vs_non*annotated.tsv "$i"_vs_non*ORFs.tsv | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1 | sort -u > "$i"_step1_speciespecific_nonORFan.txt
grep -vf "$i"_step1_speciespecific_nonORFan.txt "$i"_step1_genusspecific_ORFan.txt > "$i"_step1_speciespecific_ORFan.txt
faSomeRecords "$i"_protein_queryfile.faa "$i"_step1_speciespecific_ORFan.txt "$i"_step1_speciesspecific_ORFan.faa
done

for i in Ecoli Salmonella Mycobacterium
do
cut -f1 -d "(" "$i"_step1_genusspecific_ORFan.txt | sed "s/.*/\"&\"/g" | grep -F -f - "$i"_queryfile.gtf | gtf2bed | bedtools getfasta -s -name -fi ../Ecoli_list/"$i"_all_genomes.faa -bed - | sed "s/BALROG/balrog/g" | sed "s/GMS2/gms2/g" | sed "s/PRODIGAL/prodigal/g" | sed "s/SMORFER/smorfer/g" > "$i"_step1_genusspecific_ORFan.CDS.faa
done

#PROTEIN SEARCH PHASE - DONE#

#FLANK PHASE

for i in Ecoli Salmonella Mycobacterium
do
grep -i "prodigal" ../Ecoli_list/"$i"_annotated.gtf | bedtools sort -i | sed "s/PRODIGAL/prodigal/g" > "$i"_prodigal_annotated.gtf
done

#Get they flanks
#proximal
for i in Ecoli Salmonella Mycobacterium
do
cut -f1 -d "(" "$i"_step1_genusspecific_ORFan.txt | sed "s/.*/\"&\"/g" | grep -F -f - "$i"_queryfile.gtf | bedtools sort -i - | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$4,$6,$7,$8,$9}' |  awk -F '\t' '($4>0)' | sed "s/ \"/ \"left500_/g" | gtf2bed | bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/"$i"_all_genomes.faa -bed - > "$i"_proxflanks.faa
cut -f1 -d "(" "$i"_step1_genusspecific_ORFan.txt | sed "s/.*/\"&\"/g" | grep -F -f - "$i"_queryfile.gtf | bedtools sort -i - | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$5+500,$6,$7,$8,$9}' | sed "s/ \"/ \"right500_/g" | gtf2bed | bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/"$i"_all_genomes.faa -bed - >> "$i"_proxflanks.faa
done

#geneflanks
for j in Ecoli Salmonella Mycobacterium
do
for i in $(cut -f1 -d "(" "$j"_step1_genusspecific_ORFan.txt)
do
echo $i | sed "s/.*/\"&\"/g" | grep -F -f - "$j"_queryfile.gtf | bedtools closest -a - -b "$j"_prodigal_annotated.gtf -D a -iu -io | awk -F '\t' '($12=="CDS")' | cut -f 10-18 | sed "s/ \"/ \""$i"_down_/g" | sed 's/""/"/g' | sed 's/"_/_/g' | gtf2bed | bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/"$j"_all_genomes.faa -bed - >> "$j"_geneflanks.faa
echo $i | sed "s/.*/\"&\"/g" | grep -F -f - "$j"_queryfile.gtf | bedtools closest -a - -b "$j"_prodigal_annotated.gtf -D a -id -io | awk -F '\t' '($12=="CDS")' | cut -f 10-18 | sed "s/ \"/ \""$i"_up_/g" | sed 's/""/"/g' | sed 's/"_/_/g' | gtf2bed | bedtools getfasta -s -name -fi /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/"$j"_all_genomes.faa -bed - >> "$j"_geneflanks.faa
done
done

#Blastn proxflanks and geneflanks
for i in prox gene
do
blastn -query Ecoli_"$i"flanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Ecoli_"$i"flanks_extragenus_blastn
blastn -query Ecoli_"$i"flanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Ecoli_"$i"flanks_intragenus_blastn
blastn -query Salmonella_"$i"flanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Salmonella_"$i"flanks_extragenus_blastn
blastn -query Salmonella_"$i"flanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Enterica_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Salmonella_"$i"flanks_intragenus_blastn
blastn -query Mycobacterium_"$i"flanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Mycobacterium_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Mycobacterium_"$i"flanks_extragenus_blastn
blastn -query Mycobacterium_"$i"flanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Tuberculosis_excluded.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Mycobacterium_"$i"flanks_intragenus_blastn
blastn -query Ecoli_"$i"flanks.faa -db /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Ecoli_"$i"flanks_pangenome_blastn
blastn -query Salmonella_"$i"flanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_pangenome.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Salmonella_"$i"flanks_pangenome_blastn
blastn -query Mycobacterium_"$i"flanks.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Mycobacterium_pangenome.genomes -outfmt '6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps sstrand qstart qend sstart send qlen slen evalue bitscore' -num_threads 104 -max_target_seqs 100000 -max_hsps 1 -out Mycobacterium_"$i"flanks_pangenome_blastn
done

#targets:
for i in Ecoli Salmonella Mycobacterium
do
cat "$i"_proxflanks_extragenus_blastn "$i"_proxflanks_intragenus_blastn "$i"_proxflanks_pangenome_blastn | cut -f-2 | sed "s/_/\t/" | sed "s/\t/%/2" | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | egrep "left.*right|right.*left" | cut -f1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > "$i"_proxflanks_targets.txt
cat "$i"_geneflanks_extragenus_blastn "$i"_geneflanks_intragenus_blastn "$i"_geneflanks_pangenome_blastn | sed "s/_up_/%up%/g" | sed "s/_down_/%down%/g" | cut -f -2 | sed "s/\t/%/" | awk -F '%' '{print $2"\t"$1"%"$4}' | sort -u | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | egrep "up.*down|down.*up" | cut -f1 | awk -F'%' '{ values[$1] = (values[$1] == "" ? $2 : values[$1] ", " $2) } END { for (value in values) { print value "\t" values[value] } }' > "$i"_geneflanks_targets.txt
done

#targetlist and intervalinfo:

#Ecoli
for i in $(cut -f 1 Ecoli_proxflanks_targets.txt -d "(")
do
echo $i | sed "s/$/(/g" | grep -f - Ecoli_proxflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > Ecoli_"$i"_proxflanks_targetlist.txt
grep -w -F -f Ecoli_"$i"_proxflanks_targetlist.txt Ecoli_prox*blastn | cut -f2- -d ":" | grep "$i(" | cut -f2- -d "_" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > Ecoli_"$i"_proxflanks_intervalinfo
done
#Salmonella
for i in $(cut -f 1 Salmonella_proxflanks_targets.txt -d "(")
do
echo $i | sed "s/$/(/g" | grep -f - Salmonella_proxflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > Salmonella_"$i"_proxflanks_targetlist.txt
grep -w -F -f Salmonella_"$i"_proxflanks_targetlist.txt Salmonella_prox*blastn | cut -f2- -d ":" | grep "$i(" | cut -f2- -d "_" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > Salmonella_"$i"_proxflanks_intervalinfo
done
#Mycobacterium
for i in $(cut -f 1 Mycobacterium_proxflanks_targets.txt -d "(")
do
echo $i | sed "s/$/(/g" | grep -f - Mycobacterium_proxflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > Mycobacterium_"$i"_proxflanks_targetlist.txt
grep -w -F -f Mycobacterium_"$i"_proxflanks_targetlist.txt Mycobacterium_prox*blastn | cut -f2- -d ":" | grep "$i(" | cut -f2- -d "_" | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > Mycobacterium_"$i"_proxflanks_intervalinfo
done

#geneflanks:
for i in Ecoli Salmonella Mycobacterium
do
cat "$i"_gene*blastn | sort -u > "$i"_geneflanks_allblastn #for later
done

#Ecoli
for i in $(cut -f 1 Ecoli_geneflanks_targets.txt)
do
echo $i | sed "s/$/\t/g" | grep -F -f - Ecoli_geneflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > Ecoli_"$i"_geneflanks_targetlist.txt
cat Ecoli_geneflanks_allblastn | grep -w -F -f Ecoli_"$i"_geneflanks_targetlist.txt - | grep "$i"_ | sed "s/_up_/\t/g" | sed "s/_down_/\t/g" | cut -f1,3- | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > Ecoli_"$i"_geneflanks_intervalinfo
done
#Salmonella
for i in $(cut -f 1 Salmonella_geneflanks_targets.txt)
do
#echo $i | sed "s/$/\t/g" | grep -F -f - Salmonella_geneflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 | sort -k1 > Salmonella_"$i"_geneflanks_targetlist.txt
sort -k2 Salmonella_geneflanks_allblastn | join -1 2 -2 1 - Salmonella_"$i"_geneflanks_targetlist.txt | sed "s/ /\t/g" | awk -F '\t' '{OFS=FS}{print $2,$0}' | cut -f 1,2,4- | grep "$i"_ | sed "s/_up_/\t/g" | sed "s/_down_/\t/g" | cut -f1,3- | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > Salmonella_"$i"_geneflanks_intervalinfo
done
#Mycobacterium
for i in $(cut -f 1 Mycobacterium_geneflanks_targets.txt)
do
#echo $i | sed "s/$/\t/g" | grep -F -f - Mycobacterium_geneflanks_targets.txt | sed "s/,/\n/g" | sed "s/\t/\n/g" | sed "s/^ *//g" | tail -n+2 > Mycobacterium_"$i"_geneflanks_targetlist.txt
cat Mycobacterium_geneflanks_allblastn | grep -w -F -f Mycobacterium_"$i"_geneflanks_targetlist.txt - | grep "$i"_ | sed "s/_up_/\t/g" | sed "s/_down_/\t/g" | cut -f1,3- | awk -F '\t' '{OFS=""}{print $13,"%",$14,"%",$10,"\t",$2}' | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/%plus, /%/g" | sed "s/%minus, /%/g" | sed "s/%plus/\tplus/g" | sed "s/%minus/\tminus/g" | sed "s/%/,/g" | sed "s/\t/,/g" | sed "s/ //g" | awk -F',' '{identifier = $1; values = $2 "," $3 "," $4 "," $5; split(values, array, ","); asort(array); middle1 = array[2]; middle2 = array[3]; difference = middle2 - middle1; if (difference >= 0) { print identifier, middle1, middle2, difference, $6; } else { print identifier, middle2, middle1, -difference, $6; } }' > Mycobacterium_"$i"_geneflanks_intervalinfo
done

for j in Ecoli Salmonella Mycobacterium
do
for i in $(cut -f1 -d "(" ../"$j"_step1_speciespecific_ORFan.txt)
do
cat "$j"_"$i"_geneflanks_intervalinfo "$j"_"$i"_proxflanks_intervalinfo | grep -v -w -F -f /stor/scratch/Ochman/hassan/100724_Complete_Genomes/"$j"_pangenome_ids.txt - > temp
value=$(echo $i | sed "s/.*/\"&\"/g" | grep -F -f - ../"$j"_queryfile.gtf | awk -F '\t' '{print $5-$4+1}')
sed "s/$/ $value/" temp | awk '(($4>($6*0.5))&&($4<10000))' > "$j"_"$i"_compiled.intervalinfo.txt
done
done

for j in Ecoli Salmonella Mycobacterium
do
for i in $(grep -v -w -F -f ../"$j"_step1_speciespecific_ORFan.txt ../"$j"_step1_genusspecific_ORFan.txt | cut -f1 -d "(")
do
cat "$j"_"$i"_geneflanks_intervalinfo "$j"_"$i"_proxflanks_intervalinfo | grep -v -w -F -f /stor/scratch/Ochman/hassan/100724_Complete_Genomes/"$j"_pangenome_ids.txt - | grep -v -w -F -f /stor/scratch/Ochman/hassan/100724_Complete_Genomes/"$j"_intragenus_ids.txt - > temp
value=$(echo $i | sed "s/.*/\"&\"/g" | grep -F -f - ../"$j"_queryfile.gtf | awk -F '\t' '{print $5-$4+1}')
sed "s/$/ $value/" temp | awk '(($4>($6*0.5))&&($4<10000))' > "$j"_"$i"_compiled.intervalinfo.txt
done
done

for variable in $(ls *_compiled.intervalinfo.txt | rev | cut -f 2- -d "_" | rev)
do
    command="grep \"plus\" ${variable}_compiled.intervalinfo.txt | awk '{OFS=\"\"}{print \"samtools faidx /stor/scratch/Ochman/hassan/100724_Complete_Genomes/individual_genomes/\",\$1,\" \",\$1,\":\",\$2,\"-\",\$3}' | sed \"s/\$/ >> ${variable}_interval.faa/g\" && grep \"minus\" ${variable}_compiled.intervalinfo.txt | awk '{OFS=\"\"}{print \"samtools faidx -i /stor/scratch/Ochman/hassan/100724_Complete_Genomes/individual_genomes/\",\$1,\" \",\$1,\":\",\$2,\"-\",\$3}' | sed \"s/\$/ >> ${variable}_interval.faa/g\""
    output_file="${variable}_joint_samtools.sh"
    bash -c "$command" > "$output_file"
done

for i in $(ls *_interval.faa | rev | cut -f2- -d "_" | rev)
do
makeblastdb -in "$i"_interval.faa -dbtype nucl -out "$i"_interval
done

for j in Ecoli Salmonella Mycobacterium
do
for i in $(ls "$j"*_interval.faa | cut -f2- -d "_" | rev | cut -f2- -d "_" | rev)
do
ls "$j"_"$i"_interval.faa | cut -f2- -d "_" | rev | cut -f2- -d "_" | rev | sed "s/$/(/g" | grep --no-group-separator -A1 -f - ../"$j"_step1_genusspecific_ORFan.CDS.faa |
blastn -query - -db "$j"_"$i"_interval -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out "$j"_"$i"_interval_blastn -word_size 7
done
done

export PERL5LIB=/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/lib/
for i in $(ls *_interval_blastn | rev | cut -f3- -d "_" | rev)
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_interval_blastn > ${i}_blastn_mviewed"
done > running.sh

for j in Ecoli Salmonella Mycobacterium
do
for i in $(ls "$j"*_blastn_mviewed | cut -f2- -d "_" | rev | cut -f3- -d "_" | rev)
do
querylength=$(echo $i | sed "s/$/\(/g" | grep -A1 -f - ../"$j"_step1_genusspecific_ORFan.CDS.faa | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 "$j"_"$i"_blastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+9 "$j"_"$i"_blastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>50)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | sed "s/>/>flank_/g" > "$j"_"$i"_blastn_seq.faa
tail -n+8 "$j"_"$i"_blastn_mviewed | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 >> "$j"_"$i"_blastn_seq.faa
fi
done
done












#length of total thing
#get rid of pangenome
#get rid of species




for i in Ecoli Salmonella Mycobacterium
do
mkdir "$i"_speciesspecific_flanks
mkdir "$i"_genusspecific_flanks
cut -f1 -d "(" ../"$i"_step1_speciespecific_ORFan.txt | sed "s/^/"$i"_/g" | sed "s/$/_*intervalinfo "$i"_speciesspecific_flanks/g" | sed "s/^/cp /g" | bash
grep -v -F -f ../"$i"_step1_speciespecific_ORFan.txt ../"$i"_step1_genusspecific_ORFan.txt | cut -f1 -d "(" | sed "s/^/"$i"_/g" | sed "s/$/_*intervalinfo "$i"_genusspecific_flanks/g" | sed "s/^/cp /g" | bash
done

cd Ecoli_speciespecific_flanks
for i in $(ls *intervalinfo | rev | cut -f3- -d "_" | rev | sort -u); do ls "$i"_*intervalinfo | sed "s/^/cat /g" | bash >> "$i"_compiled_intervalinfo.txt; done
#Get rid of pangenome IDs:
for i in $(ls *_compiled_intervalinfo.txt | rev | cut -f3- -d "_" | rev | sed "s/Ecoli_//g")
do
grep "@" ../Ecoli_intervalinfo_taxonomy.tsv | cut -f1 | grep -v -w -F -f - Ecoli_"$i"_compiled_intervalinfo.txt > Ecoli_"$i"_compiled_intervalinfo.pangenomeremoved.txt
done

find . -type f -empty -delete

#Remove very short or very long interval sequences
for i in $(ls *_compiled_intervalinfo.pangenomeremoved.txt | cut -f2- -d "_" | rev | cut -f3- -d "_" | rev)
do
value=$(echo $i | sed "s/Ecoli_//g" | sed "s/.*/\"&\"/g" | grep -F -f - ../../Ecoli_queryfile.gtf | awk -F '\t' '{print $5-$4+1}')
sed "s/$/ $value/" Ecoli_"$i"_compiled_intervalinfo.pangenomeremoved.txt | awk '(($4>($6*0.5))&&($4<10000))' > Ecoli_"$i"_compiled_intervalinfo.final.txt
done

#Extract sequences:
for variable in $(ls Ecoli*_compiled_intervalinfo.final.txt | cut -f2- -d "_" | rev | cut -f 3- -d "_" | rev)
do
    command="grep \"plus\" Ecoli_${variable}_compiled_intervalinfo.final.txt | awk '{OFS=\"\"}{print \"samtools faidx /stor/scratch/Ochman/hassan/100724_Complete_Genomes/individual_genomes/\",\$1,\" \",\$1,\":\",\$2,\"-\",\$3}' | sed \"s/\$/ >> Ecoli_${variable}_interval.faa/g\" && grep \"minus\" Ecoli_${variable}_compiled_intervalinfo.final.txt | awk '{OFS=\"\"}{print \"samtools faidx -i /stor/scratch/Ochman/hassan/100724_Complete_Genomes/individual_genomes/\",\$1,\" \",\$1,\":\",\$2,\"-\",\$3}' | sed \"s/\$/ >> Ecoli_${variable}_interval.faa/g\""
    output_file=Ecoli_"${variable}_joint_samtools.sh"
    bash -c "$command" > "$output_file"
done

conda deactivate
ls *_joint_samtools.sh | sed "s/^/bash /g" > running.sh

for i in $(ls *_interval.faa | cut -f2- -d "_" | rev | cut -f2- -d "_" | rev)
do
makeblastdb -in Ecoli_"$i"_interval.faa -dbtype nucl -out Ecoli_"$i"_interval
done

for i in $(ls *_interval.faa | cut -f2- -d "_" | rev | cut -f2- -d "_" | rev)
do
ls Ecoli_"$i"_interval.faa | cut -f2- -d "_" | rev | cut -f2- -d "_" | rev | sed "s/$/(/g" | grep --no-group-separator -A1 -f - ../../Ecoli_step1_genusspecific_ORFan.CDS.faa |
blastn -query - -db Ecoli_"$i"_interval -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Ecoli_"$i"_interval_blastn -word_size 7
done

export PERL5LIB=/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/lib/
for i in $(ls *_interval_blastn | rev | cut -f3- -d "_" | rev)
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_interval_blastn > ${i}_blastn_mviewed"
done > running.sh

for j in Ecoli Salmonella Mycobacterium
do
for i in $(ls "$j"*_blastn_mviewed | cut -f2- -d "_" | rev | cut -f3- -d "_" | rev)
do
querylength=$(echo $i | sed "s/$/\(/g" | grep -A1 -f - ../"$j"_step1_genusspecific_ORFan.CDS.faa | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 "$j"_"$i"_blastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+9 "$j"_"$i"_blastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>50)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear > "$j"_"$i"_blastn_seq.faa
tail -n+8 "$j"_"$i"_blastn_mviewed | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 >> "$j"_"$i"_blastn_seq.faa
fi
done



#############REGULAR BLASTN################
#Regular blastn:

blastn -query Ecoli_step1_genusspecific_ORFan.CDS.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Escherichia_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Ecoli_extragenus_regular_blastn
blastn -query Ecoli_step1_genusspecific_ORFan.CDS.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Escherichia_db/Ecoli_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Ecoli_intragenus_regular_blastn
blastn -query Salmonella_step1_genusspecific_ORFan.CDS.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Salmonella_extragenus_regular_blastn
blastn -query Salmonella_step1_genusspecific_ORFan.CDS.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Enterica_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Salmonella_intragenus_regular_blastn
blastn -query Mycobacterium_step1_genusspecific_ORFan.CDS.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Mycobacterium_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Mycobacterium_extragenus_regular_blastn
blastn -query Mycobacterium_step1_genusspecific_ORFan.CDS.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Tuberculosis_excluded.genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Mycobacterium_intragenus_regular_blastn
blastn -query Ecoli_step1_genusspecific_ORFan.CDS.faa -db /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Ecoli_pangenome_regular_blastn
blastn -query Salmonella_step1_genusspecific_ORFan.CDS.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_pangenome.genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Salmonella_pangenome_regular_blastn
blastn -query Mycobacterium_step1_genusspecific_ORFan.CDS.faa -db /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/Mycobacterium_pangenome.genomes -outfmt 0 -num_threads 72 -num_descriptions 1000000 -num_alignments 1000000 -evalue 200000 -out Mycobacterium_pangenome_regular_blastn

mkdir regular_blastn
cp ../*genus*regular* .

for i in Ecoli Salmonella Mycobacterium
do
for j in intra extra
do
awk -v prefix="$i"_"$j"_ '/^Query=/ {close(file); match($0, / (.*)\(/, arr); file = prefix arr[1] "_blastn"} {if (file) print > file}' "$i"_"$j"genus_regular_blastn
done
done

grep "No hits found" *_blastn | grep -v "regular" | cut -f1 -d ":" | sed "s/^/rm /g" | bash

#Ecoli, extra:

for i in Ecoli Salmonella Mycobacterium
do
for j in intra extra
do
header=$(head -14 "$i"_"$j"genus_regular_blastn)
footer=$(tail -11 "$i"_"$j"genus_regular_blastn)
for file in $(ls "$i"*"$j"_*blastn); do
    # Prepend the header to the file
    { echo "$header"; cat "$file"; } > temp_file && mv temp_file "$file"
    
    # Append the footer to the file
    { cat "$file"; echo "$footer"; } > temp_file && mv temp_file "$file"
done
done
done

export PERL5LIB=/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/lib/
for i in $(ls *_blastn | grep -v "regular" | rev | cut -f2- -d "_" | rev)
do
echo "/stor/scratch/Ochman/hassan/genomics_toolbox/mview-1.67/bin/mview -in blast ${i}_blastn > ${i}_blastn_mviewed"
done > running.sh

cat ../Ecoli_CDS_queryfile.faa ../Salmonella_CDS_queryfile.faa ../Mycobacterium_CDS_queryfile.faa | seqkit fx2tab | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > all_CDS_queryfile.faa

for i in $(ls *mviewed | grep -v "regular" | rev | cut -f3- -d "_" | rev)
do
querylength=$(echo $i | cut -f3- -d "_" | sed "s/$/(/g" | grep -F -A1 -f - all_CDS_queryfile.faa | grep -v "^>" | awk '{print length($0)}')
ratio=$(tail -n+8 "$i"_blastn_mviewed | head -1 | awk '{print $(NF-1)}' | sed "s/:/\t/g" | awk -v var=$querylength -F '\t' '{print ($2-$1+1)/var}')
if (( $(echo "$ratio > 0.5" | bc -l) ))
then
tail -n+9 "$i"_blastn_mviewed | head -n-3 | awk '{print $2,$(NF-4),$NF}' | sed "s/%//g" | awk '($2>50)' | awk '{print $1,$3}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear > "$i"_blastn_seq.faa
tail -n+8 "$i"_blastn_mviewed | head -n-3 | sed "1s/^/g /g" | awk '{print $2,$NF}' | sed "s/^/>/g" | sed "s/ /\n/g" | linear | head -2 | seqkit fx2tab | cut -f2- -d "@" | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" >> "$i"_blastn_seq.faa
fi
done

#Filter by genus and species-specific sequences

#For species specific ones:
for j in Ecoli Salmonella Mycobacterium
do
for i in $(cut -f1 -d "(" ../"$j"_step1_speciespecific_ORFan.txt)
do
head -n-2 "$j"_extra_"$i"_blastn_seq.faa | sed "s/>/>regular_/g" > "$j"_"$i"_blastn_seq.regular.compiled.faa
head -n-2 "$j"_intra_"$i"_blastn_seq.faa | sed "s/>/>regular_/g" >> "$j"_"$i"_blastn_seq.regular.compiled.faa
cat "$j"_intra_"$i"_blastn_seq.faa "$j"_extra_"$i"_blastn_seq.faa | tail -2 >> "$j"_"$i"_blastn_seq.regular.compiled.faa
done
done

#For genus specific ones only:
for j in Ecoli Salmonella Mycobacterium
do
for i in $(grep -v -F -f ../"$j"_step1_speciespecific_ORFan.txt ../"$j"_step1_genusspecific_ORFan.txt | cut -f1 -d "(")
do
head -n-2 "$j"_extra_"$i"_blastn_seq.faa | sed "s/>/>regular_/g" > "$j"_"$i"_blastn_seq.regular.compiled.faa
cat "$j"_intra_"$i"_blastn_seq.faa "$j"_extra_"$i"_blastn_seq.faa | tail -2 >> "$j"_"$i"_blastn_seq.regular.compiled.faa
done
done

#Putting flanks and regulars together:
mkdir /stor/work/Ochman/hassan/protogene_extension/comparative_genomics/all_noncoding_tracing
cd /stor/work/Ochman/hassan/protogene_extension/comparative_genomics/all_noncoding_tracing
ls *_blastn_seq.faa /stor/work/Ochman/hassan/protogene_extension/comparative_genomics/regular_blastn/*_blastn_seq.regular.compiled.faa | sed "s/\/stor\/work\/Ochman\/hassan\/protogene_extension\/comparative_genomics\/regular_blastn\///g" | sed "s/_blastn_seq.faa//g" | sed "s/_blastn_seq.regular.compiled.faa//g" | sort -u | awk '{print "/stor/work/Ochman/hassan/protogene_extension/comparative_genomics/flanks/"$0"_blastn_seq.faa","/stor/work/Ochman/hassan/protogene_extension/comparative_genomics/regular_blastn/"$0"_blastn_seq.regular.compiled.faa >> "$0"_noncoding.interim"}' | sed "s/^/cat /g" | bash

for j in Ecoli Salmonella Mycobacterium
do
for i in $(ls "$j"_*noncoding.interim | cut -f2- -d "_" | rev | cut -f2- -d "_" | rev)
do
seqkit fx2tab "$j"_"$i"_noncoding.interim | sed "s/(+)//g" | sed "s/(-)//g" | sort -u | tail -n+2 | sed "s/_/\t/" | sed "s/\t$//g" | sed "s/:/\t/" | awk -F'\t' '$0 ~ /^regular/ {$2 = $2 FS "FILLER"} 1' | sed "s/ /\t/g" | sort -k2 | join -1 2 -2 1 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/"$j"_all_contig_taxonomy.tsv | sed "s/ /\t/g" | awk -F'\t' '!seen[$5,$4]++' | awk -F '\t' '{print $5":"$1":"$3":"$2"\t"$4}' | sed "s/^/>/g" | sed "s/\t/\n/g" > "$j"_"$i"_noncoding.nr
tail -2 "$j"_"$i"_noncoding.interim >> "$j"_"$i"_noncoding.nr
done
done

for i in $(ls *_noncoding.nr | rev | cut -f2- -d "_" | rev)
do
mafft --auto "$i"_noncoding.nr > "$i"_mafft.aln
done
















#DENOVO - LATER

#1. Collapse all varieties of prox and gene flanks into one file per gene cluster:

for i in $(ls *intervalinfo | rev | cut -f3- -d "_" | rev | sort -u); do ls "$i"_*intervalinfo | sed "s/^/cat /g" | bash >> "$i"_compiled_intervalinfo.txt; done


#2. Add in names to each intervalinfo file using the file all_contig_protein_taxonomy.tsv, which has been prepared using the code in assigning_conservation_to_genes.sh
cd /stor/work/Ochman/hassan/protogene_extension/comparative_genomics/flanks
cat Ecoli_*info | cut -f1 -d " " | sort -u | sort -k1 | join -1 1 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_contig_taxonomy.intervalinfo.tsv | sort -u | sed "s/ /\t/g" | awk -F '\t' '($2!="Escherichia")' > Ecoli_intervalinfo_taxonomy.tsv

for i in $(ls *_compiled_intervalinfo.txt | rev | cut -f3- -d "_" | rev | sort -u)
do
sort -k1 "$i"_compiled_intervalinfo.txt -o "$i"_compiled_intervalinfo.txt
cut -f1 -d " " "$i"_compiled_intervalinfo.txt | sort -u > temp
grep -w -F -f temp Ecoli_intervalinfo_taxonomy.tsv | sort -k1 | join -1 1 -2 1 - "$i"_compiled_intervalinfo.txt | sed 's/ [^ ]*@/ Ecoli@/' > "$i"_compiled_intervalinfo.taxa.txt
done

#3. Tagging with taxonomic information/presence-absence:

#Intra-genus:
cat Ecoli_vs_noncoliEscherichia_annotated.tsv | awk -F '\t' '($5>60&&$16<0.001)' | grep -F -f Ecoli_step1_genusspecific_ORFan.txt | cut -f-2 | rev | cut -f2- -d "_" | rev > Ecoli_intragenus_distribution.interim.txt
cat Ecoli_vs_noncoliEscherichia_ORFs.tsv | awk -F '\t' '($5>60&&$16<0.001)' | grep -F -f Ecoli_step1_genusspecific_ORFan.txt | cut -f -2 | rev | cut -f2- -d "." | rev >> Ecoli_intragenus_distribution.interim.txt
#Pangenome:
cat Ecoli_vs_pangenome_annotated.tsv | awk -F '\t' '($5>60&&$16<0.001)' | grep -F -f Ecoli_step1_genusspecific_ORFan.txt | cut -f1 -d "@" > Ecoli_pangenome_distribution.interim.txt
cat Ecoli_vs_pangenome_ORFs.tsv | awk -F '\t' '($5>60&&$16<0.001)' | grep -F -f Ecoli_step1_genusspecific_ORFan.txt | rev | cut -f2- -d "_" | rev >> Ecoli_pangenome_distribution.interim.txt

#Put them together:
cat Ecoli_intragenus_distribution.interim.txt Ecoli_pangenome_distribution.interim.txt | sort -k2 | join -1 2 -2 2 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Ecoli_intragenus_pangenome_contig_taxa.tsv | sort -u | sed 's/ [^ ]*@/ Ecoli@/' | sed "s/ /\t/g" > Ecoli_intragenus_pangenone_presence_absence.tsv

for i in $(ls *_compiled_intervalinfo.txt | rev | cut -f3- -d "_" | rev | sort -u)
do
value=$(echo $i | sed "s/Ecoli_//g" | sed "s/.*/\"&\"/g" | grep -F -f - ../Ecoli_queryfile.gtf | awk -F '\t' '{print $5-$4+1}')
sed -i "s/$/ $value/" "$i"_compiled_intervalinfo.taxa.txt
echo $i | sed "s/Ecoli_//g" | sed "s/$/\(/g" | grep -f - ../Ecoli_intragenus_pangenone_presence_absence.tsv | cut -f3 | grep -v -w -F -f - "$i"_compiled_intervalinfo.taxa.txt | awk '(($5>($7*0.5))&&($5<10000))' > "$i"_compiled_intervalinfo.taxa.final.txt
done

find . -type f -empty -delete

#Figure out different routes of analysis for different gene sets

1. E-coli proto-genes: de novo analysis

2. All other categories: Traceability to outgroups (outside genus for genus-specific, outside species for species-specific)
  a) Ecoli proto-genes
  b) Ecoli, Myc, Sal species-specific genes
  c) Ecoli, Myc, Sal genus-specific genes

#Let's do Ecoli denovo genes first to get one thing out of the way
#Tag mafft input files with these:

#Genus-specific



