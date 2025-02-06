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
Ecoli:
#Salmonella:
#cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/
#/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/all_proteins.faa --db Salmonella_pangenome_proteins
#cat /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/genomes/*_genomic.fna > /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/Salmonella_pangenome.genomes.faa
#getorf -sequence /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/Salmonella_pangenome.genomes.faa -outseq /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_db/Salmonella_pangenome.genomes.faa -table 1 -minsize 30 -find 3
seqkit fx2tab Escherichia_excluded.genomes.getorf.all | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/" > Escherichia_excluded.genomes.getorf.ATG_TTG_GTG
/stor/work/Ochman/hassan/tools/faTrans -stop Escherichia_excluded.genomes.getorf.ATG_TTG_GTG Escherichia_excluded.genomes.getorf.ATG_TTG_GTG.prot.faa

#Mycobacterium:
cd /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_db/
/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond makedb --in /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Mycobacterium_queries/Mycobacterium_331/all_331_proteins.faa --db Mycobacterium_pangenome_proteins
cat /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/genomes/*_genomic.fna > /stor/scratch/Ochman/hassan/100724_Complete_Genomes/Salmonella_queries/Salmonella_pangenome.genomes.faa

#Across pangenome, annotated
#/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_annotated.tsv -k 0 -b8 -c1
#Across pangenome, ORFs
#/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_ORFs.tsv -k 0 -b8 -c1

for i in Ecoli Salmonella Mycobacterium
do
cat "$i"_vs_GBRS_annotated.tsv "$i"_vs_GBRS_ORFs.tsv | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1 | sort -u > "$i"_step1_genusspecific_nonORFan.txt
grep -vf "$i"_step1_genusspecific_nonORFan.txt "$i"_protein_queryfile.faa | grep "^>" | tr -d ">" > "$i"_step1_genusspecific_ORFan.txt
faSomeRecords "$i"_protein_queryfile.faa "$i"_step1_genusspecific_ORFan.txt "$i"_step1_genusspecific_ORFan.faa
done

#Get they flanks

