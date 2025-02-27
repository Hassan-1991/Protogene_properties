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

mkdir 042323_ostir_calculation
cat x*out | awk -F ',' '($3==51)' | grep "^WP" | sed "s/$/,annotated/g" > ostir_plotting.tsv
cat x*out | awk -F ',' '($3==51)' | grep "^sORF" | sed "s/$/,Riboseq/g" >> ostir_plotting.tsv
cat x*out | awk -F ',' '($3==51)' | grep "^REL606" | grep -v "other" | sed "s/$/,MS/g" >> ostir_plotting.tsv
cat x*out | awk -F ',' '($3==51)' | grep "^REL606" | grep "other" | sed "s/$/,others/g" >> ostir_plotting.tsv

sed -i "1s/^/name,start_codon,start_position,expression,RBS_distance_bp,dG_total,dG_rRNA:mRNA,dG_mRNA,dG_spacing,dG_standby,dG_start_codon,type\n/g" ostir_plotting.tsv
