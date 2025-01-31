#Mori2021

#RNAseq
#Those without Rifampicin treatment
#accession list: downloaded from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA847230&o=acc_s%3Aa
for i in $(cat accession_list.txt); do esearch -db sra -query $i | efetch -format runinfo | cut -f1,12 -d ','; done | grep -v "LibraryName" > SRR_GSM.csv
#GSM_filename.interim: copy-pasted from here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205717
#filename_conditions.interim: copy-paste from supp file 3. science.abk2066_tables_s3_to_s7
sed "s/\t/,/g" filename_conditions.interim | sort -k1 -t ',' > filename_conditions.csv
sort -k2 GSM_filename.interim | sed "s/\t/,/g" | join -t ',' -1 2 -2 1 - filename_conditions.csv | sort -t ',' -k2 > GSM_filename_conditions.csv
sort -k2 -t ',' SRR_GSM.csv | join -t ',' -1 2 -2 2 - GSM_filename_conditions.csv > SRR_GSM_filename_conditions.csv

#MS

#Caglar 2017

#RNAseq

#MS
