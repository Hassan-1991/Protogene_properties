#This code generates a transcription-supported protein database for two MS datasets

#Caglar2017 dataset - REL606
#Mori2021 - K12 MG1655, K12 NCM3722. One is from Nissle, which we can ignore

#Are these strains sufficiently similar to use ORFs generated from just one database?

cp /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/sequence_RS.fasta K12MG1655.faa
cp /stor/work/Ochman/hassan/proteomics_denovo/04092024_compgenomics_final/GCF_000017985.1_ASM1798v1_genomic.fna REL606.faa

#Caglar2017, MS mgf files:
#Each MURI_XXX directory contains all corresponding mgf files
/stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/
#Caglar2017, RNAseq bam files:
/stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/RNAseq/bamfiles/
#Mori2021, MS raw files:
cd /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/
#Get the index file from pride depository first:
wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD014948/
#Get file names, append to ftp directory name, add wget
grep -P "wiff" *html | cut -f8 -d "\"" | sort -u | sort -u | grep -v "SW" | sed "s/^/https:\/\/ftp.pride.ebi.ac.uk\/pride\/data\/archive\/2021\/03\/PXD014948\//g" | sed "s/^/wget /g" | split -l 10 -
#parallelize, 10 per screen:
ls x* | sed "s/^/bash /g" > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh
#Mori2021, RNAseq raw files:
cd /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/RNAseq/raw
#Download accession_list.txt from here: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA847230&o=acc_s%3Aa&s=SRR19588406,SRR19588407,SRR19588408,SRR19588409,SRR19588410,SRR19588411,SRR19588412,SRR19588413,SRR19588414,SRR19588415,SRR19588416,SRR19588417,SRR19588418,SRR19588419,SRR19588420,SRR19588421,SRR19588422,SRR19588423,SRR19588424,SRR19588425,SRR19588426,SRR19588427,SRR19588428,SRR19588429,SRR19588430,SRR19588431,SRR19588432,SRR19588433,SRR19588434,SRR19588435,SRR19588436,SRR19588437,SRR19588438,SRR19588439,SRR19588440,SRR19588441,SRR19588442,SRR19588443,SRR19588444,SRR19588445,SRR19588446,SRR19588447,SRR19588448,SRR19588449,SRR19588450
#Fetch SRR files:
sed "s/^/\/stor\/work\/Ochman\/hassan\/tools\/sratoolkit.3.0.6-ubuntu64\/bin\/prefetch /g" accession_list.txt | bash
#Get fastq out of them:
ls SRR*/*sra | sed "s/^/\/stor\/work\/Ochman\/hassan\/tools\/sratoolkit\.3\.0\.6-ubuntu64\/bin\/fasterq-dump /g" | bash
#Mori2021, RNAseq bam files:
#QC:
ls *fastq | sed "s/^/fastp -i /g" | sed "s/$/ -o /g" | awk '{print $0$3"trimmed"}' | sed "s/fastqtrimmed/trimmed.fastq/g" | sed "s/$/ --thread 16/g" | bash
##
