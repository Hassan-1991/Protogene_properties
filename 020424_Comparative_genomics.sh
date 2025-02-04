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

#Across pangenome, annotated
#/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_proteins --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_annotated.tsv -k 0 -b8 -c1
#Across pangenome, ORFs
#/stor/work/Ochman/hassan/E.coli_ORFan/E.coli_ORFan_pipeline_8-10/diamond blastp -q Ecoli_protein_queryfile.faa -d /stor/work/Ochman/hassan/Ecoli_pangenome/103024_updated_pipeline/backup/all_450_genomes.getorf.ATG_TTG_GTG.prot --outfmt 6 qseqid sseqid pident nident qcovhsp length mismatch gapopen gaps qstart qend sstart send qlen slen evalue bitscore --ultra-sensitive --out all_proteins_vs_pangenome_ORFs.tsv -k 0 -b8 -c1
