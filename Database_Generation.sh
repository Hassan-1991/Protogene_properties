#This code generates proteogenomics search databases from genomes

#Caglar2017 dataset was generated from REL606
#All but one experiment in Mori2021 was conducted in NCM3722, which is 99%+ similar to K12_MG1655; so we just use that

cd /stor/scratch/Ochman/hassan/112724_protogene_extension/proteogenomics_database

#Place strain genomes in current directory
#These genomes were downloaded as a part of protogene_curation

cp /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/sequence_RS.fasta NC_000913.3.faa
cp /stor/work/Ochman/hassan/proteomics_denovo/04092024_compgenomics_final/GCF_000017985.1_ASM1798v1_genomic.fna NC_012967.1.faa
/stor/work/Ochman/hassan/proteomics_denovo/8925_raw_data/ECOR_11_genome.faa .
/stor/work/Ochman/hassan/proteomics_denovo/8925_raw_data/ECOR_27_genome.faa .
/stor/work/Ochman/hassan/proteomics_denovo/8925_raw_data/ECOR_37_genome.faa .

#Annotate using prodigal, balrog, genemarks2
#For some reason balrog doesn't run for the ECOR genomes
for i in NC_000913.3.faa NC_012967.1.faa ECOR_11_genome.faa ECOR_27_genome.faa ECOR_37_genome.faa
do echo "/stor/work/Ochman/hassan/tools/prod_gms_balrog.sh ${i}"
done > running.sh

#annotate with smorfinder

conda activate smorfinder_env
for i in NC_000913.3.faa NC_012967.1.faa ECOR_11_genome.faa ECOR_27_genome.faa ECOR_37_genome.faa
do
smorf single "$i" && mv smorf_output "$i"_smorf_output
done

#convert to gtf

for i in NC_000913.3.faa NC_012967.1.faa ECOR_11_genome.faa ECOR_27_genome.faa ECOR_37_genome.faa
do
cut -f-8 "$i"*smorf_output/*gff | awk -F '\t' '{print $0,"\ttranscript_id \""$1"_smorfer_"NR"\";gene_id \""$1"_smorfer_"NR"\";"}' > "$i"_smorfer.gtf
cut -f-8 "$i"*prodigal.gff | awk -F '\t' '{print $0,"\ttranscript_id \""$1"_prodigal_"NR"\";gene_id \""$1"_prodigal_"NR"\";"}' > "$i"_prodigal.gtf
awk -F '\t' '($3=="CDS")' "$i"*genemarks2.gff | cut -f -8 | awk -F '\t' '{print $0,"\ttranscript_id \""$1"_gms2_"NR"\";gene_id \""$1"_gms2_"NR"\";"}' > "$i"_gms2.gtf
done

#Concatenate everything

for i in NC_000913.3 NC_012967.1 ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
cat "$i"*smorfer.gtf "$i"*prodigal.gtf "$i"*gms2.gtf "$i"*balrog.gtf  | grep -v "^#" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,".",$7,".",$9,$10}' > "$i"_annotated.gtf
done

for i in NC_000913.3 NC_012967.1 ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
getorf -sequence "$i".faa -outseq "$i".getorf.bacterial -table 11 -minsize 30 -find 3
done

for i in NC_000913.3 NC_012967.1 ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
grep "^>" "$i".getorf.bacterial | grep -v "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$2,"\t",$4+3,"\t.\t+\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- > "$i".getorf.bacterial.gtf
grep "^>" "$i".getorf.bacterial | grep "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$4-3,"\t",$2,"\t.\t-\t0\tgene_id \"",$1,"\";transcript_id \"",$1,"\";"}' | sed "s/_/\t/2" | cut -f1,3- >> "$i".getorf.bacterial.gtf
done

cp NC_012967.1.faa REL606.faa
cp NC_012967.1_annotated.gtf REL606_annotated.gtf
cp NC_012967.1.getorf.bacterial.gtf REL606.getorf.bacterial.gtf

cp NC_000913.3.faa K12MG1655.faa
cp NC_000913.3_annotated.gtf K12MG1655_annotated.gtf
cp NC_000913.3.getorf.bacterial.gtf K12MG1655.getorf.bacterial.gtf

for i in K12MG1655 REL606 ECOR_11_genome ECOR_27_genome ECOR_37_genome
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

for i in K12MG1655 REL606 ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
sed "s/(+)//g" "$i".getorf.bacterial.longestORFlengths.tsv | sed "s/(-)//g" | sort -k1 > temp
cut -f -2 -d "\"" "$i".getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="+")' | awk '{print $2"\t.\tCDS\t"$6-$NF+1"\t"$6"\t.\t+\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$5}' | awk -F '\t' '($5-$4+1>32)' > "$i".getorf.bacterial.longest.gtf
cut -f -2 -d "\"" "$i".getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="-")' | awk '{print $2"\t.\tCDS\t"$5"\t"$5+$NF-1"\t.\t-\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$6}' | awk -F '\t' '($5-$4+1>32)' >> "$i".getorf.bacterial.longest.gtf
awk -F '\t' '($7=="+")' "$i"_annotated.gtf | cut -f5 | sort -u | sed "s/^/\t/g" | sed "s/$/\t./g" | grep -o -F -f - "$i".getorf.bacterial.longest.gtf > "$i"_stops_represented
awk -F '\t' '($7=="-")' "$i"_annotated.gtf | cut -f4 | sort -u | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | grep -o -F -f - "$i".getorf.bacterial.longest.gtf >> "$i"_stops_represented
grep -v -F -f "$i"_stops_represented "$i"_annotated.gtf | awk -F '\t' '($7=="+")' | cut -f5 | sed "s/^/\t/g" | sed "s/$/\t./g" | sort -u > "$i"_unrepresented_stops
grep -v -F -f "$i"_stops_represented "$i"_annotated.gtf | awk -F '\t' '($7=="-")' | cut -f4 | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | sort -u >> "$i"_unrepresented_stops
cut -f 1-9 "$i".getorf.bacterial.longest.gtf | bedtools intersect -s -wo -a - -b "$i"_annotated.gtf | awk -F '\t' '(($7=="+")&&($5==$14))' | cut -f2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > "$i"_exclude
cut -f 1-9 "$i".getorf.bacterial.longest.gtf | bedtools intersect -s -wo -a - -b "$i"_annotated.gtf | awk -F '\t' '(($7=="-")&&($4==$13))' | cut -f2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" >> "$i"_exclude
cut -f 1-9 "$i".getorf.bacterial.longest.gtf | grep -v -F -f "$i"_exclude > "$i".getorf.bacterial.longest.annotexcluded.gtf
cat "$i".getorf.bacterial.longest.annotexcluded.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > "$i".getorf.bacterial.longest.final.faa
cat "$i".getorf.bacterial.longest.annotexcluded.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - "$i".getorf.bacterial.longest.annotexcluded.gtf > "$i".getorf.bacterial.longest.final.gtf
#Non-redundantify the annotated gtf
awk -F '\t' '($7=="+")' "$i"_annotated.gtf | cut -f5 | sort -u > "$i"_annotated_stops_plus
for j in $(cat "$i"_annotated_stops_plus); do awk -v var="$j" -F '\t' '($7=="+"&&$5==var)' "$i"_annotated.gtf | sort -nk 4 | head -1; done > "$i"_annotated.nr.gtf
awk -F '\t' '($7=="-")' "$i"_annotated.gtf | cut -f4 | sort -u > "$i"_annotated_stops_minus
for j in $(cat "$i"_annotated_stops_minus); do awk -v var="$j" -F '\t' '($7=="-"&&$4==var)' "$i"_annotated.gtf | sort -nrk 5 | head -1; done >> "$i"_annotated.nr.gtf
cat "$i"_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > "$i"_all_annotated.nr.seqkit.faa
cat "$i"_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - "$i"_annotated.nr.gtf > "$i"_all_annotated.nr.seqkit.gtf
cat "$i"_all_annotated.nr.seqkit.gtf "$i".getorf.bacterial.longest.final.gtf | sort -nk4 > "$i".final.gtf
cat "$i".final.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - > "$i".final.faa
/stor/work/Ochman/hassan/tools/faTrans -stop "$i".final.faa "$i".final.prot.faa
done
