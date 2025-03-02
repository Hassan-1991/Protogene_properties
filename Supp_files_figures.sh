#Fig S1	Annotated proteins detected across different conditions, inc. for those w just one peptide
#Table S1	List of all proteins detected by MS
#Table S2	list of studies w inclusion exclusion criteria
#Table S3	Proteins of different categories. Maybe the ORFan + traceability + intergenics table
#Fig S2	ORFan gene bar graph for both categories
#Fig S3	Taxonomic conservation of conserved Novel vs annotated genes
#Fig S4	Expression pattern broken down between ORFan and conserved proto-genes
#Fig S5	Annotated gene pangenome
#Table S4	Pangenome distribution for all species-specific genes
#Fig S6	Consistency bar plot between datasets

#Let's check consistency of proteins across the two datasets
#Start by mapping REL606 and K12 proteins together

######

/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d K12MG1655 -f 2 --gff3-fasta-annotation=1 ../REL606.final.faa > REL606_mappedto_K12MG1655.gff3
awk -F '\t' '($3=="mRNA")' REL606_mappedto_K12MG1655.gff3 | grep "mrna1" | cut -f1 -d ';' | sed "s/ID=//g" | sed "s/.mrna1$//g" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' > REL606_mappedto_K12MG1655.gtf

egrep -i "prodigal|gms2|balrog|smorf" K12MG1655.final.gtf > K12.annot.gtf
egrep -i "prodigal|gms2|balrog|smorf" REL606_mappedto_K12MG1655.gtf > REL606.annot.mappedtoK12.gtf

awk '
    $3 == "CDS" {
        stop = ($7 == "+") ? $5 : $4;
        key = $1 ":" stop ":" $7;
        gene = gensub(/.*gene_name "([^"]+)".*/, "\\1", "g", $0);
        genes[key] = (key in genes) ? genes[key] "," gene : gene;
    }
    END {
        for (k in genes)
            if (genes[k] ~ /,/) print k, genes[k];
    }
' K12.annot.gtf REL606.annot.mappedtoK12.gtf | sort -k1,1 > K12.gtf.interim

awk -F '\t' '($7=="+")' K12.gtf.interim | cut -f5,9,17 | sed "s/\t/\"/1" | cut -f1,3,7 -d "\"" | sed "s/\"/\t/g" > K12.REL606.interim2
awk -F '\t' '($7=="-")' K12.gtf.interim | cut -f4,9,17 | sed "s/\t/\"/1" | cut -f1,3,7 -d "\"" | sed "s/\"/\t/g" >> K12.REL606.interim2

#Now let's count up the number of times a protein showed up

awk -F '\t' '($16<0.0001)' ../data/Caglar2017/MS_spectra_searches/MURI*tsv | egrep -i "prodigal|smorf|gms2|balrog" | cut -f 10,11 | sort -u | grep -v "XXX" > annotated_proteins_detected.tsv
awk -F '\t' '($16<0.0001)' ../RNAseq_diversion/MURI*tsv | egrep -i "prodigal|smorf|gms2|balrog" | cut -f 10,11 | sort -u | grep -v "XXX" >> annotated_proteins_detected.tsv
awk -F '\t' '($16<0.0001)' ../data/Mori2021/MS/mgf/real_mgf/chludwig*tsv | egrep -i "prodigal|smorf|gms2|balrog" | cut -f 10,11 | sort -u | grep -v "XXX" >> annotated_proteins_detected.tsv

awk 'NR==FNR {map[$3] = $2; next} {if ($2 in map) $2 = map[$2]; print}' K12.REL606.interim2 annotated_proteins_detected.tsv | grep -v "REL606" | sed "s/ /\t/g" | sed "s/(+)//g" | sed "s/(-)//g" > annotated_protein_peptides_qval0001.tsv

cut -f1 annotated_protein_peptides_qval0001.tsv | cut -f2- -d "." | rev | cut -f2- -d "." | rev | sed "s/^/./g" | sed "s/$/./g" | grep -F -f - ../data/Caglar2017/MS_spectra_searches/MURI*tsv | awk -F '\t' '($16<0.01)' | grep -v "XXX" > annotated_results.interim
cut -f1 annotated_protein_peptides_qval0001.tsv | cut -f2- -d "." | rev | cut -f2- -d "." | rev | sed "s/^/./g" | sed "s/$/./g" | grep -F -f - ../RNAseq_diversion/MURI*tsv | awk -F '\t' '($16<0.01)' | grep -v "XXX" >> annotated_results.interim
cut -f1 annotated_protein_peptides_qval0001.tsv | cut -f2- -d "." | rev | cut -f2- -d "." | rev | sed "s/^/./g" | sed "s/$/./g" | grep -F -f - ../data/Mori2021/MS/mgf/real_mgf/chludwig*tsv | awk -F '\t' '($16<0.01)' | grep -v "XXX" >> annotated_results.interim

cut -f2 annotated_proteins_detected.tsv | grep -F -f - annotated_results.interim > annotated_results.interim.2 #get the peptides
awk -F '\t' 'NR==FNR {map[$3] = $2; next} {if ($11 in map) $11 = map[$11]; print}' OFS='\t' K12.REL606.interim2 annotated_results.interim.2 | grep -v "REL606" > annotated_results.fixed.tsv #replace protein names
sed "s/(+)//g" annotated_results.fixed.tsv | sed "s/(-)//g" > temp && mv temp annotated_results.fixed.tsv

cat annotated_results.fixed.tsv | sed "s/(+)//g" | sed "s/(-)//g" | awk -F '\t' '($16<0.01)' | cut -f1,11 | rev | cut -f1 -d ":" | rev | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.01"}' > annot_protein_distribution.csv
cat annotated_results.fixed.tsv | sed "s/(+)//g" | sed "s/(-)//g" | awk -F '\t' '($16<0.001)' | cut -f1,11 | rev | cut -f1 -d ":" | rev | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.001"}' >> annot_protein_distribution.csv
cat annotated_results.fixed.tsv | sed "s/(+)//g" | sed "s/(-)//g" | awk -F '\t' '($16<0.0001)' | cut -f1,11 | rev | cut -f1 -d ":" | rev | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.0001"}' >> annot_protein_distribution.csv
sort -t ',' -k1 annot_protein_distribution.csv -o annot_protein_distribution.csv

cat annotated_results.fixed.tsv | sed "s/(+)//g" | sed "s/(-)//g" | awk -F '\t' '($16<0.01)' | cut -f1,11 | rev | cut -f1 -d ":" | rev | sort -u | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.01"}' > annot_SAMPLENO_distribution.csv
cat annotated_results.fixed.tsv | sed "s/(+)//g" | sed "s/(-)//g" | awk -F '\t' '($16<0.01)' | cut -f1,11 | rev | cut -f1 -d ":" | rev | sort -u | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.01"}' >> annot_SAMPLENO_distribution.csv
cat annotated_results.fixed.tsv | sed "s/(+)//g" | sed "s/(-)//g" | awk -F '\t' '($16<0.01)' | cut -f1,11 | rev | cut -f1 -d ":" | rev | sort -u | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.01"}' >> annot_SAMPLENO_distribution.csv
sort -t ',' -k1 annot_SAMPLENO_distribution.csv -o annot_SAMPLENO_distribution.csv

#TO GET MEDIAN VALUUES FROM BOTH:
awk '{a[NR]=$1} END {print (NR%2 ? a[(NR+1)/2] : (a[NR/2] + a[NR/2+1]) / 2)}'

cut -f2 annotated_protein_peptides_qval0001.tsv | sort | uniq -c | grep " 1 " | rev | cut -f1 -d " " | rev | sort -u > just1pept.txt

######
#all spectra:
#Extract list of all protein and their spectra identification + amino acid sequence + file

 /stor/scratch/Ochman/hassan/112724_protogene_extension/promising_spectra/*spectra.tsv
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/Caglar2017_RNAseq_promising_spectra/*spectra.tsv
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/promising_spectra/*spectra.tsv
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR2023_promising_spectra/*spectra.tsv

tail -n+2 protein_PSM_distribution.csv | cut -f2 -d ',' | cut -f1 -d "(" | sort -u > 39_proteins.txt

#Get protein/peptide/sample/scan info:
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/promising_spectra/*spectra.tsv \
   /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/Caglar2017_RNAseq_promising_spectra/*spectra.tsv \
   /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/promising_spectra/*spectra.tsv \
   /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR2023_promising_spectra/*spectra.tsv |
grep -F -f 39_proteins.txt - | while read -r file; do
    # Extract parts from the filename
    part1=$(basename "$file" | rev | cut -f1 -d "/" | rev | cut -f1-3 -d "_")
    part2=$(basename "$file" | rev | cut -f1 -d "/" | rev | cut -f4- -d "_" | rev | cut -f3- -d "." | rev)

    # Extract data from the file itself while preserving spaces
    title=$(grep "TITLE=" "$file" | sed "s/TITLE=//g" | paste -sd ",")
    seq=$(grep "SEQ=" "$file" | sed "s/SEQ=//g" | paste -sd ",")

    # Use printf to ensure proper tab separation while keeping spaces intact
    printf "%s\t%s\t%s\t%s\n" "$part1" "$part2" "$title" "$seq"
done



#Random. Let's count the number of rbs in two genomes
#We can start with longest_bacterial

for i in Ecoli Salmonella Mycobacterium
do
#getorf -sequence "$i".focalstrain.faa -outseq "$i".getorf.bacterial -table 1 -minsize 30 -find 3
seqkit fx2tab "$i".getorf.bacterial | sed "s/\t$//g" | grep -P -iv "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/g" | grep "^>" | grep -v "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$2,"\t",$4+3,"\t.\t+\t0\ttranscript_id \"",$1,"\";gene_id \"",$1,"\";"}' | sed "s/_/\t/1" | cut -f1,3- > rbs_promoters/"$i".getorf.bacterial.gtf
seqkit fx2tab "$i".getorf.bacterial | sed "s/\t$//g" | grep -P -iv "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/g" | grep "^>" | grep "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$4-3,"\t",$2,"\t.\t-\t0\ttranscript_id \"",$1,"\";gene_id \"",$1,"\";"}' | sed "s/_/\t/1" | cut -f1,3- >> rbs_promoters/"$i".getorf.bacterial.gtf
done

sed -i "s/^NC/NC_999999.9/g" Mycobacterium.getorf.bacterial.gtf
sed -i "s/^NC/NC_000913.3/g" Ecoli.getorf.bacterial.gtf

cat Ecoli.getorf.bacterial.gtf | awk -F '\t' '($7=="+")' | awk -F '\t' '{OFS=FS}{print "NC_000913.3",$2,$3,$4-50,$4+50,$6,$7,$8,$9}' | awk '($4>0)' | gtf2bed | bedtools getfasta -s -name -fi ../Ecoli.focalstrain.faa -bed - > Ecoli_ostir_input.faa
cat Ecoli.getorf.bacterial.gtf | awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print "NC_000913.3",$2,$3,$5-50,$5+50,$6,$7,$8,$9}' | awk '($4>0)' | gtf2bed | bedtools getfasta -s -name -fi ../Ecoli.focalstrain.faa -bed - >> Ecoli_ostir_input.faa

cat Mycobacterium.getorf.bacterial.gtf | awk -F '\t' '($7=="+")' | awk -F '\t' '{OFS=FS}{print "NC_999999.9",$2,$3,$4-50,$4+50,$6,$7,$8,$9}' | awk '($4>0)' | gtf2bed | bedtools getfasta -s -name -fi ../Mycobacterium.focalstrain.faa -bed - > Mycobacterium_ostir_input.faa
cat Mycobacterium.getorf.bacterial.gtf | awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print "NC_999999.9",$2,$3,$5-50,$5+50,$6,$7,$8,$9}' | awk '($4>0)' | gtf2bed | bedtools getfasta -s -name -fi ../Mycobacterium.focalstrain.faa -bed - >> Mycobacterium_ostir_input.faa

cat Ecoli_ostir_input.faa | split -l10000 - Ecoli_
cat Mycobacterium_ostir_input.faa | split -l10000 - Mycobacterium_

ls Ecoli_a* Mycobacterium_a* | awk '{print "ostir -i "$0" -o "$0".ostir -j 104"}' > running.sh
###RUN FROM HERE###

cat Ecoli*ostir | awk -F ',' '($3==51)' > Ecoli.ostir
cat Mycobacterium*ostir | awk -F ',' '($3==51)' > Mycobacterium.ostir

#Doesn't seem like Myc has more than EC, so ignore.
