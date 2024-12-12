#This code generates proteogenomics search databases from two genomes

#Caglar2017 dataset was generated from REL606
#All but one experiment in Mori2021 was conducted in NCM3722, which is 99%+ similar to K12_MG1655; so we just use that

cd /stor/scratch/Ochman/hassan/112724_protogene_extension/proteogenomics_database

#Place strain genomes in current directory
#These genomes were downloaded as a part of protogene_curation

cp /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/sequence_RS.fasta NC_000913.3.faa
cp /stor/work/Ochman/hassan/proteomics_denovo/04092024_compgenomics_final/GCF_000017985.1_ASM1798v1_genomic.fna NC_012967.1.faa

#Extract strain-specific annotation files
#Ecoli_all_annotated.gtf was generated using the protogene_curation code

for i in NC_012967.1 NC_000913.3
do 
sed "s/REL606/NC_012967.1/g" /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/Ecoli_all_annotated.gtf | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,$6=".",$7,$8,$9}' | awk -F '\t' '($3=="CDS")' | grep -F "$i" > "$i".all_annotated.gtf
done

for i in NC_012967.1 NC_000913.3
do
getorf -sequence "$i".faa -outseq "$i".getorf.bacterial -table 11 -minsize 30 -find 3
done



grep "^>" REL606.getorf.bacterial | grep -v "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$2,"\t",$4+3,"\t.\t+\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- > REL606.getorf.bacterial.gtf
grep "^>" REL606.getorf.bacterial | grep "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$4-3,"\t",$2,"\t.\t-\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- >> REL606.getorf.bacterial.gtf

grep "^>" K12MG1655.getorf.bacterial | grep -v "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$2,"\t",$4+3,"\t.\t+\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- > K12MG1655.getorf.bacterial.gtf
grep "^>" K12MG1655.getorf.bacterial | grep "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$4-3,"\t",$2,"\t.\t-\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- >> K12MG1655.getorf.bacterial.gtf

awk -F '\t' '($7=="+")' REL606.getorf.bacterial.gtf |
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$5,$6,$7,$8,$9}' |
awk -F '\t' '{OFS=FS}{if ($4 < 1) $4 = 1; print}' |
gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - |
sed '/^>/!{ s/.*/echo "&" | rev/e }' |
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' |
cut -f2- -d " " |
sed -E 's/AGT|AAT|GAT/%/g' |
seqkit fx2tab | cut -f1 -d "%" |
sed -E 's/GTA|GTG|GTT/@/g' | grep "@" |
rev | cut -f2- -d "@" | rev |
sed "s/@/NNN/g" | sed "s/ //g" | awk -F '\t' '{print $1,length($2)+6}' | sed "s/ /\t/g" > REL606.getorf.bacterial.longestORFlengths.tsv

awk -F '\t' '($7=="-"&&$4>0)' REL606.getorf.bacterial.gtf |
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+500,$6,$7,$8,$9}' |
awk -F '\t' '{OFS=FS}{if ($5 > 4629812) $5 = 4629812; print}' |
gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - |
sed '/^>/!{ s/.*/echo "&" | rev/e }' |
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' |
cut -f2- -d " " |
sed -E 's/AGT|AAT|GAT/%/g' |
seqkit fx2tab | cut -f1 -d "%" |
sed -E 's/GTA|GTG|GTT/@/g' | grep "@" |
rev | cut -f2- -d "@" | rev |
sed "s/@/NNN/g" | sed "s/ //g" | awk -F '\t' '{print $1,length($2)+6}' | sed "s/ /\t/g" >> REL606.getorf.bacterial.longestORFlengths.tsv

awk -F '\t' '($7=="+")' K12MG1655.getorf.bacterial.gtf |
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-500,$5,$6,$7,$8,$9}' |
awk -F '\t' '{OFS=FS}{if ($4 < 1) $4 = 1; print}' |
gtf2bed | bedtools getfasta -s -name -fi K12MG1655.faa -bed - |
sed '/^>/!{ s/.*/echo "&" | rev/e }' |
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' |
cut -f2- -d " " |
sed -E 's/AGT|AAT|GAT/%/g' |
seqkit fx2tab | cut -f1 -d "%" |
sed -E 's/GTA|GTG|GTT/@/g' | grep "@" |
rev | cut -f2- -d "@" | rev |
sed "s/@/NNN/g" | sed "s/ //g" | awk -F '\t' '{print $1,length($2)+6}' | sed "s/ /\t/g" > K12MG1655.getorf.bacterial.longestORFlengths.tsv

awk -F '\t' '($7=="-"&&$4>0)' K12MG1655.getorf.bacterial.gtf |
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+500,$6,$7,$8,$9}' |
awk -F '\t' '{OFS=FS}{if ($5 > 4641652) $5 = 4641652; print}' |
gtf2bed | bedtools getfasta -s -name -fi K12MG1655.faa -bed - |
sed '/^>/!{ s/.*/echo "&" | rev/e }' |
awk '{if ($0 ~ /^>/) print; else {gsub(/.{3}/, "& "); print}}' |
cut -f2- -d " " |
sed -E 's/AGT|AAT|GAT/%/g' |
seqkit fx2tab | cut -f1 -d "%" |
sed -E 's/GTA|GTG|GTT/@/g' | grep "@" |
rev | cut -f2- -d "@" | rev |
sed "s/@/NNN/g" | sed "s/ //g" | awk -F '\t' '{print $1,length($2)+6}' | sed "s/ /\t/g" >> K12MG1655.getorf.bacterial.longestORFlengths.tsv

sed "s/(+)//g" REL606.getorf.bacterial.longestORFlengths.tsv | sed "s/(-)//g" | sort -k1 > temp
sed "s/(+)//g" K12MG1655.getorf.bacterial.longestORFlengths.tsv | sed "s/(-)//g" | sort -k1 > temp

cut -f -2 -d "\"" REL606.getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="+")' | awk '{print $2"\t.\tCDS\t"$6-$NF+1"\t"$6"\t.\t+\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$5}' | awk -F '\t' '($5-$4+1>32)' > REL606.getorf.bacterial.longest.gtf
cut -f -2 -d "\"" REL606.getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="-")' | awk '{print $2"\t.\tCDS\t"$5"\t"$5+$NF-1"\t.\t-\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$6}' | awk -F '\t' '($5-$4+1>32)' >> REL606.getorf.bacterial.longest.gtf

cut -f -2 -d "\"" K12MG1655.getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="+")' | awk '{print $2"\t.\tCDS\t"$6-$NF+1"\t"$6"\t.\t+\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$5}' | awk -F '\t' '($5-$4+1>32)' > K12MG1655.getorf.bacterial.longest.gtf
cut -f -2 -d "\"" K12MG1655.getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="-")' | awk '{print $2"\t.\tCDS\t"$5"\t"$5+$NF-1"\t.\t-\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$6}' | awk -F '\t' '($5-$4+1>32)' >> K12MG1655.getorf.bacterial.longest.gtf

awk -F '\t' '($7=="+")' REL606_all_annotated.mod.gtf | cut -f5 | sort -u | sed "s/^/\t/g" | sed "s/$/\t./g" | grep -o -F -f - REL606.getorf.bacterial.longest.gtf > stops_represented
awk -F '\t' '($7=="-")' REL606_all_annotated.mod.gtf | cut -f4 | sort -u | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | grep -o -F -f - REL606.getorf.bacterial.longest.gtf >> stops_represented

awk -F '\t' '($7=="+")' K12MG1655_all_annotated.mod.gtf | cut -f5 | sort -u | sed "s/^/\t/g" | sed "s/$/\t./g" | grep -o -F -f - K12MG1655.getorf.bacterial.longest.gtf > stops_represented
awk -F '\t' '($7=="-")' K12MG1655_all_annotated.mod.gtf | cut -f4 | sort -u | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | grep -o -F -f - K12MG1655.getorf.bacterial.longest.gtf >> stops_represented

grep -v -F -f stops_represented REL606_all_annotated.mod.gtf | awk -F '\t' '($7=="+")' | cut -f5 | sed "s/^/\t/g" | sed "s/$/\t./g" | sort -u > unrepresented_stops
grep -v -F -f stops_represented REL606_all_annotated.mod.gtf | awk -F '\t' '($7=="-")' | cut -f4 | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | sort -u >> unrepresented_stops

grep -v -F -f stops_represented K12MG1655_all_annotated.mod.gtf | awk -F '\t' '($7=="+")' | cut -f5 | sed "s/^/\t/g" | sed "s/$/\t./g" | sort -u > unrepresented_stops
grep -v -F -f stops_represented K12MG1655_all_annotated.mod.gtf | awk -F '\t' '($7=="-")' | cut -f4 | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | sort -u >> unrepresented_stops

cut -f 1-9 REL606.getorf.bacterial.longest.gtf | bedtools intersect -s -wo -a - -b REL606_all_annotated.mod.gtf | awk -F '\t' '(($7=="+")&&($5==$14))' | cut -f2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > exclude
cut -f 1-9 REL606.getorf.bacterial.longest.gtf | bedtools intersect -s -wo -a - -b REL606_all_annotated.mod.gtf | awk -F '\t' '(($7=="-")&&($4==$13))' | cut -f2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" >> exclude

cut -f 1-9 K12MG1655.getorf.bacterial.longest.gtf | bedtools intersect -s -wo -a - -b K12MG1655_all_annotated.mod.gtf | awk -F '\t' '(($7=="+")&&($5==$14))' | cut -f2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > exclude
cut -f 1-9 K12MG1655.getorf.bacterial.longest.gtf | bedtools intersect -s -wo -a - -b K12MG1655_all_annotated.mod.gtf | awk -F '\t' '(($7=="-")&&($4==$13))' | cut -f2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" >> exclude

cut -f 1-9 REL606.getorf.bacterial.longest.gtf | grep -v -F -f exclude > REL606.getorf.bacterial.longest.annotexcluded.gtf
cat REL606.getorf.bacterial.longest.annotexcluded.gtf | gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > REL606.getorf.bacterial.longest.final.faa
cat REL606.getorf.bacterial.longest.annotexcluded.gtf | gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - REL606.getorf.bacterial.longest.annotexcluded.gtf > REL606.getorf.bacterial.longest.final.gtf

cut -f 1-9 K12MG1655.getorf.bacterial.longest.gtf | grep -v -F -f exclude > K12MG1655.getorf.bacterial.longest.annotexcluded.gtf
cat K12MG1655.getorf.bacterial.longest.annotexcluded.gtf | gtf2bed | bedtools getfasta -s -name -fi K12MG1655.faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > K12MG1655.getorf.bacterial.longest.final.faa
cat K12MG1655.getorf.bacterial.longest.annotexcluded.gtf | gtf2bed | bedtools getfasta -s -name -fi K12MG1655.faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - K12MG1655.getorf.bacterial.longest.annotexcluded.gtf > K12MG1655.getorf.bacterial.longest.final.gtf

#Non-redundantify the annotated gtf

awk -F '\t' '($7=="+")' REL606_all_annotated.mod.gtf | cut -f5 | sort -u > annotated_stops_plus
for i in $(cat annotated_stops_plus); do awk -v var="$i" -F '\t' '($7=="+"&&$5==var)' REL606_all_annotated.mod.gtf | sort -nk 4 | head -1; done > REL606_all_annotated.nr.gtf
awk -F '\t' '($7=="-")' REL606_all_annotated.mod.gtf | cut -f4 | sort -u > annotated_stops_minus
for i in $(cat annotated_stops_minus); do awk -v var="$i" -F '\t' '($7=="-"&&$4==var)' REL606_all_annotated.mod.gtf | sort -nrk 5 | head -1; done >> REL606_all_annotated.nr.gtf

awk -F '\t' '($7=="+")' K12MG1655_all_annotated.mod.gtf | cut -f5 | sort -u > annotated_stops_plus
for i in $(cat annotated_stops_plus); do awk -v var="$i" -F '\t' '($7=="+"&&$5==var)' K12MG1655_all_annotated.mod.gtf | sort -nk 4 | head -1; done > K12MG1655_all_annotated.nr.gtf
awk -F '\t' '($7=="-")' K12MG1655_all_annotated.mod.gtf | cut -f4 | sort -u > annotated_stops_minus
for i in $(cat annotated_stops_minus); do awk -v var="$i" -F '\t' '($7=="-"&&$4==var)' K12MG1655_all_annotated.mod.gtf | sort -nrk 5 | head -1; done >> K12MG1655_all_annotated.nr.gtf

cat REL606_all_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > REL606_all_annotated.nr.seqkit.faa
cat REL606_all_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - REL606_all_annotated.nr.gtf > REL606_all_annotated.nr.seqkit.gtf

cat K12MG1655_all_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi K12MG1655.faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > K12MG1655_all_annotated.nr.seqkit.faa
cat K12MG1655_all_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi K12MG1655.faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - K12MG1655_all_annotated.nr.gtf > K12MG1655_all_annotated.nr.seqkit.gtf

cat REL606_all_annotated.nr.seqkit.gtf REL606.getorf.bacterial.longest.final.gtf | sort -nk4 > REL606.final.gtf
cat K12MG1655_all_annotated.nr.seqkit.gtf K12MG1655.getorf.bacterial.longest.final.gtf | sort -nk4 > K12MG1655.final.gtf

cat REL606.final.gtf | gtf2bed | bedtools getfasta -s -name -fi REL606.faa -bed - > REL606.final.faa
/stor/work/Ochman/hassan/tools/faTrans -stop REL606.final.faa REL606.final.prot.faa

cat K12MG1655.final.gtf | gtf2bed | bedtools getfasta -s -name -fi K12MG1655.faa -bed - > K12MG1655.final.faa
/stor/work/Ochman/hassan/tools/faTrans -stop K12MG1655.final.faa K12MG1655.final.prot.faa

