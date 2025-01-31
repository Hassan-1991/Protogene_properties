#Annotate ECOR11, ECOR27, ECOR37 genomes

cp /stor/work/Ochman/hassan/tools/gms2_linux_64/.gmhmmp2_key .
for i in ECOR_11_genome.faa ECOR_27_genome.faa ECOR_37_genome.faa
do echo "/stor/work/Ochman/hassan/tools/prod_gms_balrog.sh ${i}"
done > running.sh

#balrog not running for some reason

smorf single ECOR_11_genome.faa && mv smorf_output ECOR_11_genome_smorf_output

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
cut -f-8 "$i"*smorf_output/*gff | awk -F '\t' '{print $0,"\ttranscript_id \""$1"_smorfer_"NR"\";gene_id \""$1"_smorfer_"NR"\";"}' > "$i"_smorfer.gtf
cut -f-8 "$i"*prodigal.gff | awk -F '\t' '{print $0,"\ttranscript_id \""$1"_prodigal_"NR"\";gene_id \""$1"_prodigal_"NR"\";"}' > "$i"_prodigal.gtf
awk -F '\t' '($3=="CDS")' "$i"*genemarks2.gff | cut -f -8 | awk -F '\t' '{print $0,"\ttranscript_id \""$1"_gms2_"NR"\";gene_id \""$1"_gms2_"NR"\";"}' > "$i"_gms2.gtf
done

#cat em
for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
cat "$i"*smorfer.gtf "$i"*prodigal.gtf "$i"*gms2.gtf | grep -v "^#" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5,".",$7,".",$9,$10}' > "$i"_annotated.gtf
done

#Extracting ORFs:

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

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
sed "s/(+)//g" "$i".getorf.bacterial.longestORFlengths.tsv | sed "s/(-)//g" | sort -k1 > temp
cut -f -2 -d "\"" "$i".getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="+")' | awk '{print $2"\t.\tCDS\t"$6-$NF+1"\t"$6"\t.\t+\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$5}' | awk -F '\t' '($5-$4+1>32)' > "$i".getorf.bacterial.longest.gtf
cut -f -2 -d "\"" "$i".getorf.bacterial.gtf | sed "s/gene_id \"//g" | sort -k9 | join -1 9 -2 1 - temp | awk '($8=="-")' | awk '{print $2"\t.\tCDS\t"$5"\t"$5+$NF-1"\t.\t-\t.\ttranscript_id \""$1"\";gene_id \""$1"\";\t"$6}' | awk -F '\t' '($5-$4+1>32)' >> "$i".getorf.bacterial.longest.gtf
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
awk -F '\t' '($7=="+")' "$i"_annotated.gtf | cut -f5 | sort -u | sed "s/^/\t/g" | sed "s/$/\t./g" | grep -o -F -f - "$i".getorf.bacterial.longest.gtf > "$i"_stops_represented
awk -F '\t' '($7=="-")' "$i"_annotated.gtf | cut -f4 | sort -u | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | grep -o -F -f - "$i".getorf.bacterial.longest.gtf >> "$i"_stops_represented
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
grep -v -F -f "$i"_stops_represented "$i"_annotated.gtf | awk -F '\t' '($7=="+")' | cut -f5 | sed "s/^/\t/g" | sed "s/$/\t./g" | sort -u > "$i"_unrepresented_stops
grep -v -F -f "$i"_stops_represented "$i"_annotated.gtf | awk -F '\t' '($7=="-")' | cut -f4 | sed "s/^/CDS\t/g" | sed "s/$/\t/g" | sort -u >> "$i"_unrepresented_stops
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
cut -f 1-9 "$i".getorf.bacterial.longest.gtf | bedtools intersect -s -wo -a - -b "$i"_annotated.gtf | awk -F '\t' '(($7=="+")&&($5==$14))' | cut -f2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" > "$i"_exclude
cut -f 1-9 "$i".getorf.bacterial.longest.gtf | bedtools intersect -s -wo -a - -b "$i"_annotated.gtf | awk -F '\t' '(($7=="-")&&($4==$13))' | cut -f2 -d "\"" | sort -u | sed "s/.*/\"&\"/g" >> "$i"_exclude
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
cut -f 1-9 "$i".getorf.bacterial.longest.gtf | grep -v -F -f "$i"_exclude > "$i".getorf.bacterial.longest.annotexcluded.gtf
cat "$i".getorf.bacterial.longest.annotexcluded.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > "$i".getorf.bacterial.longest.final.faa
cat "$i".getorf.bacterial.longest.annotexcluded.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - "$i".getorf.bacterial.longest.annotexcluded.gtf > "$i".getorf.bacterial.longest.final.gtf
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
awk -F '\t' '($7=="+")' "$i"_annotated.gtf | cut -f5 | sort -u > "$i"_annotated_stops_plus
for j in $(cat "$i"_annotated_stops_plus); do awk -v var="$j" -F '\t' '($7=="+"&&$5==var)' "$i"_annotated.gtf | sort -nk 4 | head -1; done > "$i"_all_annotated.nr.gtf
awk -F '\t' '($7=="-")' "$i"_annotated.gtf | cut -f4 | sort -u > "$i"_annotated_stops_plus
for k in $(cat "$i"_annotated_stops_plus); do awk -v var="$k" -F '\t' '($7=="-"&&$4==var)' "$i"_annotated.gtf | sort -nrk 5 | head -1; done >> "$i"_all_annotated.nr.gtf
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
cat "$i"_all_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > "$i"_all_annotated.nr.seqkit.faa
cat "$i"_all_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - "$i"_all_annotated.nr.gtf > "$i"_all_annotated.nr.seqkit.gtf
done

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
cat "$i"_all_annotated.nr.seqkit.gtf "$i".getorf.bacterial.longest.final.gtf | sort -nk4 > "$i".final.gtf
cat "$i".final.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - > "$i".final.faa
/stor/work/Ochman/hassan/tools/faTrans -stop "$i".final.faa "$i".final.prot.faa
done

#MSGF clumsy ass. Make tsv:
cd /stor/work/Ochman/hassan/Fall_2022/massspec_database/MSGFplus/ #For some reason this code only runs in the directory of the jar file
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/*mzid | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i "$1" -o "$1".tsv -showQValue 1 -showDecoy 1 -unroll 1"}' | sed "s/mzid.tsv/tsv/g" > running.sh

for i in $(ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/*mzid | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
do
echo "python3 /stor/work/Ochman/hassan/mass_spec/test_files/mgf_search_result_annotator_test1.py --format MSGF_ident --input /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/${i}.mgf --search /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/${i}.mzid --output /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/${i}_annotated.mgf"
done > running.sh

cat 8925*tsv | egrep -v "XXX" | egrep -i "prodigal|gms|balrog|smorf" > all_canonical_PSMs.tsv
cat 8925*tsv | egrep "XXX" | cut -f10 | cut -f2- -d '.' | rev | cut -f2- -d '.' | rev | sort -u | grep -o -F -f - all_canonical_PSMs.tsv | sort -u | grep -F -f - *tsv | egrep "XXX" | cut -f11 | sort -u > exclude
cat 8925*tsv | egrep -iv "XXX|gms|balrog|smorf|prodigal" | cut -f10 | cut -f2- -d '.' | rev | cut -f2- -d '.' | rev | sort -u | grep -o -F -f - all_canonical_PSMs.tsv | sort -u | grep -F -f - *tsv | egrep -iv "XXX|gms|balrog|smorf|prodigal" | cut -f11 | sort -u >> exclude
grep -v -F -f exclude 8925*tsv | awk -F '\t' '($16<0.0001)' | cut -f11 | sort -u > test2
grep -v "REL606" test2 | grep -F -f - MURI*tsv | awk -F '\t' '($16<0.00001)' | cut -f1,11 | grep -v "XXX" > promising_spectra

egrep -vi "gms|prodigal|smorf|balrog|XXX" test2 | grep -F -f - 8925*tsv | awk -F '\t' '($16<0.00001)' | cut -f1,11 | grep -v "XXX" > promising_spectra

cut -f2 promising_spectra | sort -u | cut -f1 -d "(" | while IFS= read -r i; do
    grep -F "$i(" promising_spectra | cut -f1 -d ":" | rev | cut -f2- -d "." | rev | while IFS= read -r j; do
        grep -F "$i(" "/stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/$j.tsv" | 
        awk -F '\t' '($16<0.00001)' | cut -f4 | rev | cut -f2 -d "\"" | rev | while IFS= read -r k; do
            sed -n "/$k$/,/END IONS/p" "${j}_annotated.mgf" > "${i}_${j}.spectra"
        done
    done
done

