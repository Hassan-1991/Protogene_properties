#For Caglar2017 and ECOR datasets, convert raw files to mgf using msconvert, e.g.:
docker run -it --rm -e WINEDEBUG=-all -v `pwd`:`pwd` -w `pwd` chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /stor/work/Ochman/hassan/Fall_2022/NEW_DATA/8869_HU_01.raw --mgf
#msconvert is very finicky. It needs to be run in the same folder as the raw file. Otherwise even absolute addresses don't work

#mgf files located here:
/stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/
#ECOR mgfs:
/stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/

#For Caglar2017 RNAseq:
#Run Trimmomatic with standard options to convert all fastq files to make trimmed files
#Stock code:

java -jar ./trimmomatic-0.39.jar PE -threads 36 read1.fastq read2.fastq read1_paired.fastq read1_unpaired.fastq read2_paired.fastq read2_unpaired.fastq ILLUMINACLIP:illumina_truseq_6bp_dual_barcode_PE.fasta:2:30:10:8:true MINLEN:20

#Run bowtie to convert trimmed files to bam files
#Stock code:
bowtie2-build --threads 72 -f REL606.fasta REL606.bowtie
bowtie2 -x REL606.bowtie -1 read1_paired.fastq -2 read2_paired.fastq -U read1_unpaired.fastq,read2_unpaired.fastq -S read.sam --threads 72 --local --very-sensitive-local && samtools view -bS read.sam -@ 72 > read.bam && time samtools sort -n -o read_sorted.bam read.bam -@ 72 && time samtools index read_sorted.bam -@ 72

#bamfiles are located here:
/stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/RNAseq/bamfiles/

#For Mori2021:

#Mori2021,MS:
cd /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/
#Download MS raw files
#Get the index file from pride depository first:
wget https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD014948/
#Get file names from html file, append to ftp directory name, add wget
grep -P "wiff" *html | cut -f8 -d "\"" | sort -u | sort -u | grep -v "SW" | sed "s/^/https:\/\/ftp.pride.ebi.ac.uk\/pride\/data\/archive\/2021\/03\/PXD014948\//g" | sed "s/^/wget /g" | split -l 10 -
#parallelize, 10 per screen:
ls x* | sed "s/^/bash /g" > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh
#mgf conversion happens in Windows PC using GUI msconvert
#peak picking enabled

#Mori2021,RNAseq:
cd /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/RNAseq/raw
#Download accession_list.txt from here: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA847230&o=acc_s%3Aa&s=SRR19588406,SRR19588407,SRR19588408,SRR19588409,SRR19588410,SRR19588411,SRR19588412,SRR19588413,SRR19588414,SRR19588415,SRR19588416,SRR19588417,SRR19588418,SRR19588419,SRR19588420,SRR19588421,SRR19588422,SRR19588423,SRR19588424,SRR19588425,SRR19588426,SRR19588427,SRR19588428,SRR19588429,SRR19588430,SRR19588431,SRR19588432,SRR19588433,SRR19588434,SRR19588435,SRR19588436,SRR19588437,SRR19588438,SRR19588439,SRR19588440,SRR19588441,SRR19588442,SRR19588443,SRR19588444,SRR19588445,SRR19588446,SRR19588447,SRR19588448,SRR19588449,SRR19588450
#Fetch SRR files:
sed "s/^/\/stor\/work\/Ochman\/hassan\/tools\/sratoolkit.3.0.6-ubuntu64\/bin\/prefetch /g" accession_list.txt | bash
#Get fastq out of them:
ls SRR*/*sra | sed "s/^/\/stor\/work\/Ochman\/hassan\/tools\/sratoolkit\.3\.0\.6-ubuntu64\/bin\/fasterq-dump /g" | bash
#Trim them with fastp:
cd /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/RNAseq/fastp_trimmed
ls *fastq | sed "s/^/fastp -i /g" | sed "s/$/ -o /g" | awk '{print $0$3"trimmed"}' | sed "s/fastqtrimmed/trimmed.fastq/g" | sed "s/$/ --thread 16/g" | bash
#Generate bamfiles by mapping these to K12 genome:
cd /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/RNAseq/bamfiles
bowtie2-build --threads 72 -f K12MG1655.faa K12MG1655.bowtie
for i in $(ls ../fastp_trimmed/ | cut -f1 -d ".")
do
echo "bowtie2 -x K12MG1655.bowtie -U /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/RNAseq/fastp_trimmed/${i}.trimmed.fastq -S ${i}.sam --threads 72 --local --very-sensitive-local \&\& samtools view -bS ${i}.sam -@ 72 > ${i}.bam && time samtools sort -n -o ${i}_sorted.bam ${i}.bam -@ 72"
done

#mgf files:
/stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/
#bamfiles:
/stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/RNAseq/bamfiles

