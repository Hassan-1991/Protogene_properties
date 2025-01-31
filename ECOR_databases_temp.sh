#Annotate ECOR11, ECOR27, ECOR37 genomes

cp /stor/work/Ochman/hassan/tools/gms2_linux_64/.gmhmmp2_key .
for i in ECOR_11_genome.faa ECOR_27_genome.faa ECOR_37_genome.faa
do echo "/stor/work/Ochman/hassan/tools/prod_gms_balrog.sh ${i}"
done > running.sh

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
getorf -sequence "$i".faa -outseq "$i".getorf.bacterial -table 11 -minsize 30 -find 3
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
grep "^>" "$i".getorf.bacterial | grep -v "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$2,"\t",$4+3,"\t.\t+\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- > "$i".getorf.bacterial.gtf
grep "^>" "$i".getorf.bacterial | grep "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$4-3,"\t",$2,"\t.\t-\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- >> "$i".getorf.bacterial.gtf
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
awk -F '\t' '($7=="+")' "$i".getorf.bacterial.gtf |
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$5,$6,$7,$8,$9}' |
awk -F '\t' '{OFS=FS}{if ($4 < 1) $4 = 1; print}' |
gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - |
sed '/^>/!{ s/.*/echo "&" | rev/e }' |
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' |
cut -f2- -d " " |
sed -E 's/AGT|AAT|GAT/%/g' |
seqkit fx2tab | cut -f1 -d "%" |
sed -E 's/GTA|GTG|GTT/@/g' | grep "@" |
rev | cut -f2- -d "@" | rev |
sed "s/@/NNN/g" | sed "s/ //g" | awk -F '\t' '{print $1,length($2)+6}' | sed "s/ /\t/g" > "$i".getorf.bacterial.longestORFlengths.tsv
awk -F '\t' '($7=="-"&&$4>0)' "$i".getorf.bacterial.gtf |
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+500,$6,$7,$8,$9}' |
awk -F '\t' '{OFS=FS}{if ($5 > 4629812) $5 = 4629812; print}' |
gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - |
sed '/^>/!{ s/.*/echo "&" | rev/e }' |
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' |
cut -f2- -d " " |
sed -E 's/AGT|AAT|GAT/%/g' |
seqkit fx2tab | cut -f1 -d "%" |
sed -E 's/GTA|GTG|GTT/@/g' | grep "@" |
rev | cut -f2- -d "@" | rev |
sed "s/@/NNN/g" | sed "s/ //g" | awk -F '\t' '{print $1,length($2)+6}' | sed "s/ /\t/g" >> "$i".getorf.bacterial.longestORFlengths.tsv
done

sed "s/(+)//g" REL606.getorf.bacterial.longestORFlengths.tsv | sed "s/(-)//g" | sort -k1 > temp

cut -f -2 -d "\"" REL606.getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="+")' | awk '{print $2"\t.\tCDS\t"$6-$NF+1"\t"$6"\t.\t+\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$5}' | awk -F '\t' '($5-$4+1>32)' > REL606.getorf.bacterial.longest.gtf
cut -f -2 -d "\"" REL606.getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="-")' | awk '{print $2"\t.\tCDS\t"$5"\t"$5+$NF-1"\t.\t-\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$6}' | awk -F '\t' '($5-$4+1>32)' >> REL606.getorf.bacterial.longest.gtf

awk -F '\t' '($7=="+")' REL606_all_annotated.mod.gtf | cut -f5 | sort -u | sed "s/^/\t/g" | sed "s/$/\t./g" | grep -o -F -f - REL606.getorf.bacterial.longest.gtf > stops_represented
awk -F '\t' '($7=="-")' REL606_all_annotated.mod.gtf | cut -f4 | sort -u | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | grep -o -F -f - REL606.getorf.bacterial.longest.gtf >> stops_represented

