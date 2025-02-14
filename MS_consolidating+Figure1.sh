#Map each protein in each genome

mkdir consolidating_datasets
cd consolidating_datasets
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/REL606.faa .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/K12MG1655.faa .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_11_genome.faa .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_27_genome.faa .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_37_genome.faa .

cp /stor/scratch/Ochman/hassan/112724_protogene_extension/REL606.final.faa .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/K12MG1655.final.faa .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_11_genome.final.faa .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_27_genome.final.faa .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_37_genome.final.faa .

cp /stor/scratch/Ochman/hassan/112724_protogene_extension/REL606.final.gtf .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/K12MG1655.final.gtf .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_11_genome.final.gtf .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_27_genome.final.gtf .
cp /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/ECOR_37_genome.final.gtf .

#Build indices
for i in REL606 K12MG1655 ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap_build -D . -d "$i" "$i".faa
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 REL606.final.faa > REL606_"$i"_gmapped_ORFs.gff3
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 K12MG1655.final.faa > K12MG1655_"$i"_gmapped_ORFs.gff3
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 ECOR_11_genome.final.faa > ECOR_11_"$i"_gmapped_ORFs.gff3
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 ECOR_27_genome.final.faa > ECOR_27_"$i"_gmapped_ORFs.gff3
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d "$i" -f 2 --gff3-fasta-annotation=1 ECOR_37_genome.final.faa > ECOR_37_"$i"_gmapped_ORFs.gff3
done

#Convert to genome-specific coordinates, then see if stop codons match with searched proteins
awk -F '\t' '($3=="mRNA")' REL606_K12MG1655_gmapped_ORFs.gff3 | grep "mrna1" | cut -f1 -d ';' | sed "s/ID=//g" | sed "s/.mrna1$//g" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' > REL606_K12MG1655_gmapped_ORFs.gtf

bedtools intersect -wo -a REL606_K12MG1655_gmapped_ORFs.gtf -b K12MG1655.final.gtf | awk -F '\t' '($7==$16&&$7=="+")' | awk -F '\t' '($5==$14)' | cut -f9,18 | cut -f2,6 -d "\"" | sed "s/\"/\t/g"
bedtools intersect -wo -a REL606_K12MG1655_gmapped_ORFs.gtf -b K12MG1655.final.gtf | awk -F '\t' '($7==$16&&$7=="-")' | awk -F '\t' '($5==$14)' | cut -f9,18 | cut -f2,6 -d "\"" | sed "s/\"/\t/g"

/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d K12MG1655 -f 2 --gff3-fasta-annotation=1 ../Caglar2017_MSvalidated_proteins.cds.faa | awk -F '\t' '($3=="mRNA")' | grep "mrna1" | cut -f1 -d ";" | sed "s/ID=//g" | sed "s/.mrna1$//g" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' | bedtools intersect -wo -a - -b K12MG1655.final.gtf | awk -F '\t' '($7==$16&&$7=="+")' | awk -F '\t' '($5==$14)' | cut -f2,6 -d "\"" | sed "s/\"/\t/g" > REL606_mapped_to_K12.tsv
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d K12MG1655 -f 2 --gff3-fasta-annotation=1 ../Caglar2017_MSvalidated_proteins.cds.faa | awk -F '\t' '($3=="mRNA")' | grep "mrna1" | cut -f1 -d ";" | sed "s/ID=//g" | sed "s/.mrna1$//g" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' | bedtools intersect -wo -a - -b K12MG1655.final.gtf | awk -F '\t' '($7==$16&&$7=="-")' | awk -F '\t' '($5==$14)' | cut -f2,6 -d "\"" | sed "s/\"/\t/g" >> REL606_mapped_to_K12.tsv

grep --no-group-separator -A1 -w -F -f ECOR2023_MSvalidated_proteins.txt data/ECOR_2023/*genome.final.faa | cut -f2- -d ":" | sed '/^>/!s/.*-//' > ECOR2023_MSvalidated_proteins.cds.faa

/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d K12MG1655 -f 2 --gff3-fasta-annotation=1 ../ECOR2023_MSvalidated_proteins.cds.faa | awk -F '\t' '($3=="mRNA")' | grep "mrna1" | cut -f1 -d ";" | sed "s/ID=//g" | sed "s/.mrna1$//g" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' | bedtools intersect -wo -a - -b K12MG1655.final.gtf | awk -F '\t' '($7==$16&&$7=="+")' | awk -F '\t' '($5==$14)' | cut -f2,6 -d "\"" | sed "s/\"/\t/g" > ECOR2023_mapped_to_K12.tsv
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d K12MG1655 -f 2 --gff3-fasta-annotation=1 ../ECOR2023_MSvalidated_proteins.cds.faa | awk -F '\t' '($3=="mRNA")' | grep "mrna1" | cut -f1 -d ";" | sed "s/ID=//g" | sed "s/.mrna1$//g" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' | bedtools intersect -wo -a - -b K12MG1655.final.gtf | awk -F '\t' '($7==$16&&$7=="-")' | awk -F '\t' '($5==$14)' | cut -f2,6 -d "\"" | sed "s/\"/\t/g" >> ECOR2023_mapped_to_K12.tsv

sed "s/$/(/g" ../Mori2021 | grep --no-group-separator -A1 -F -f - K12MG1655.final.faa > Mori2021_MSvalidated_proteins.cds.faa
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d REL606 -f 2 --gff3-fasta-annotation=1 Mori2021_MSvalidated_proteins.cds.faa | awk -F '\t' '($3=="mRNA")' | grep "mrna1" | cut -f1 -d ";" | sed "s/ID=//g" | sed "s/.mrna1$//g" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' | bedtools intersect -wo -a - -b REL606.final.gtf | awk -F '\t' '($7==$16&&$7=="+")' | awk -F '\t' '($5==$14)' | cut -f2,6 -d "\"" | sed "s/\"/\t/g" > K12_mapped_to_REL606.tsv
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d REL606 -f 2 --gff3-fasta-annotation=1 Mori2021_MSvalidated_proteins.cds.faa | awk -F '\t' '($3=="mRNA")' | grep "mrna1" | cut -f1 -d ";" | sed "s/ID=//g" | sed "s/.mrna1$//g" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' | bedtools intersect -wo -a - -b REL606.final.gtf | awk -F '\t' '($7==$16&&$7=="-")' | awk -F '\t' '($5==$14)' | cut -f2,6 -d "\"" | sed "s/\"/\t/g" >> K12_mapped_to_REL606.tsv

#NC_000913.3_101381 in K12 is same as NC_012967.1_102979 in REL606
#But it's not novel, blastp search tells me it's a peptide release factor

NC_012967.1_64224(+)	NC_000913.3_64789

#I have five different genomes -  REL606, K12, ECOR11, ECOR27, ECOR37

#For each gene, there needs to be five columns - evidence in mass spec data (from absent in genome to level of evidence) 

#All Caglar2017 (REL606) proteins:

cat Caglar2017_withRNAseq Caglar2017_withoutRNAseq | sort -u | sed "s/$/(/g" | grep --no-group-separator -A1 -F -f - REL606.final.faa > Caglar2017_MSvalidated_proteins.cds.faa

#Let's see if this problem can be solved quickly by blasting
awk -F '\t' '($16<0.01)' ../data/Caglar2017/MS_spectra_searches/*tsv | grep -v "XXX_" | grep -v "REL606" | cut -f11 | sort -u > REL606_bothdatasets_qval0.01.txt
awk -F '\t' '($16<0.01)' ../RNAseq_diversion/*tsv | grep -v "XXX_" | grep -v "REL606" | cut -f11 | sort -u >> REL606_bothdatasets_qval0.01.txt
sort -u REL606_bothdatasets_qval0.01.txt -o REL606_bothdatasets_qval0.01.txt
awk -F '\t' '($16<0.01)' ../data/Mori2021/MS/mgf/real_mgf/*tsv | grep -v "XXX_" | egrep -iv "gms|balrog|smorf|prodigal" | cut -f11 | sort -u > K12MG1655_qval0.01.txt
awk -F '\t' '($16<0.01)' ../data/ECOR_2023/8925*tsv | cut -f11 | sort -u | egrep -iv "XXX|gms|balrog|prodigal|smorf" > ECOR2023_qval0.01.txt

seqkit fx2tab ../REL606.final.prot.faa | sed "s/\t$//g" | grep -f REL606_bothdatasets_qval0.01.txt - | sed "s/^/>/g" | sed "s/\t/\n/g" > REL606_bothdatasets_qval0.01.prot.faa
seqkit fx2tab ../K12MG1655.final.prot.faa | sed "s/\t$//g" | grep -f K12MG1655_qval0.01.txt - | sed "s/^/>/g" | sed "s/\t/\n/g" > K12MG1655_qval0.01.prot.faa
cat ../data/ECOR_2023/ECOR_*_genome.final.prot.faa | seqkit fx2tab - | sed "s/\t$//g" | grep -f ECOR2023_qval0.01.txt - | sed "s/^/>/g" | sed "s/\t/\n/g" > ECOR2023_qval0.01.prot.faa
cat REL606_bothdatasets_qval0.01.prot.faa K12MG1655_qval0.01.prot.faa ECOR2023_qval0.01.prot.faa > all_qval0.01.prot.faa

seqkit fx2tab all_qval0.01.prot.faa | sed "s/\t$//g" | awk -F '\t' '{print $1,length($2)}' | sed "s/ /\t/g" | sort -k1 > all_qval0.01.lengths.tsv
grep -v "^#" test | sort -k1 | join -1 1 -2 1 - all_qval0.01.lengths.tsv | sed "s/ /\t/g" | awk -F '\t' '{OFS=FS}{print $0,($8-$7+1)/$13}' | awk -F '\t' '($3>80&&$NF>0.8)' | awk -F '\t' '($1!=$2)' 

grep -v "^#" test | sort -k1 | join -1 1 -2 1 - all_qval0.01.lengths.tsv | sed "s/ /\t/g" | awk -F '\t' '{OFS=FS}{print $0,($8-$7+1)/$13}' | awk -F '\t' '($3>80&&$NF>0.8)' | awk -F '\t' '($1!=$2)' | sed "s/(+)/(/" | sed "s/(-)/(/" | grep -F -f all_MS_validated_proteins.txt - > parseable.tsv

#Parse the above later to find intra-database links

grep "^>" ../Caglar2017_MSvalidated_proteins.cds.faa | cut -f1 -d "(" | tr -d ">" | cat - ../Mori2021 ../ECOR2023_MSvalidated_proteins.txt | sort -u | sed "s/$/(/g" > all_MS_validated_proteins.txt

#Excluding annotated hits:
cut -f1 -d "(" /stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/all_MS_validated_proteins.txt | grep -v -F -f all_protogenes_tobeexcluded.txt -

#Exclude:
NC_000913.3_101381(-)
NC_000913.3_116107(-)
NC_000913.3_4041(+)
NC_012967.1_102979(-)
NC_012967.1_124406(-)

cut -f1 -d "(" /stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/all_MS_validated_proteins.txt | grep -v -F -f all_protogenes_tobeexcluded.txt - | sed "s/$/(/g" | grep -F -f - ../comparative_genomics/Ecoli_vs_pangenome_annotated.tsv | awk -F '\t' '($16<0.001&&$5>60)' | cut -f1 | sort -u | cut -f1 -d "(" | sed "s/$/(/g" | grep -F -v -f - /stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/all_MS_validated_proteins.txt | grep -v -F -f all_protogenes_tobeexcluded.txt - | cut -f1 -d "(" > /stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets/all_MS_validated_proteins.final.txt

#List all peptides
grep --no-group-separator -A1 "SEQ=" ../promising_spectra/*tsv | rev | cut -f1 -d "/" | rev | awk 'NR%2{printf "%s ", $0; next}1' | sed "s/SEQ=/:/g" | sed "s/USER03=/:/g" | sed "s/ /:/g" | cut -f1,3,5 -d ":" | sed "s/:/\t/g" > all_promising_peptides.tsv
grep --no-group-separator -A1 "SEQ=" ../RNAseq_diversion/Caglar2017_RNAseq_promising_spectra/*tsv | rev | cut -f1 -d "/" | rev | awk 'NR%2{printf "%s ", $0; next}1' | sed "s/SEQ=/:/g" | sed "s/USER03=/:/g" | sed "s/ /:/g" | cut -f1,3,5 -d ":" | sed "s/:/\t/g" >> all_promising_peptides.tsv
grep --no-group-separator -A1 "SEQ=" ../data/Mori2021/MS/mgf/real_mgf/promising_spectra/*tsv | rev | cut -f1 -d "/" | rev | awk 'NR%2{printf "%s ", $0; next}1' | sed "s/SEQ=/:/g" | sed "s/USER03=/:/g" | sed "s/ /:/g" | cut -f1,3,5 -d ":" | sed "s/:/\t/g" >> all_promising_peptides.tsv
grep --no-group-separator -A1 "SEQ=" ../data/ECOR_2023/ECOR2023_promising_spectra/*tsv | rev | cut -f1 -d "/" | rev | awk 'NR%2{printf "%s ", $0; next}1' | sed "s/SEQ=/:/g" | sed "s/USER03=/:/g" | sed "s/ /:/g" | cut -f1,3,5 -d ":" | sed "s/:/\t/g" >> all_promising_peptides.tsv

sed "s/$/_/g" all_MS_validated_proteins.final.txt | grep -F -f - all_promising_peptides.tsv | rev | cut -f1 -d "=" | rev | sort -u | grep -v "Oxidation" | cut -f2 | sed "s/C/C+57.021/g" | sort -u > 40_peptides.txt
sed "s/$/_/g" all_MS_validated_proteins.final.txt | grep -F -f - all_promising_peptides.tsv | rev | cut -f1 -d "=" | rev | sort -u | grep "Oxidation" | cut -f2- | sed "s/C/C+57.021/g" | sort -u > interim
#manually edit the interim file
cut -f1 interim | sort -u >> 40_peptides.txt

#tsv searches:
grep -F -f 40_peptides.txt /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_spectra_searches/MURI*tsv | awk -F '\t' '($16<0.01)' > all_40peptide_results.tsv
grep -F -f 40_peptides.txt /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/MURI*tsv | awk -F '\t' '($16<0.01)' >> all_40peptide_results.tsv
grep -F -f 40_peptides.txt /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/chludwig*tsv | awk -F '\t' '($16<0.01)' >> all_40peptide_results.tsv
grep -F -f 40_peptides.txt /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/8925*tsv | awk -F '\t' '($16<0.01)' >> all_40peptide_results.tsv

for i in $(cat 40_peptides.txt); do echo -n $i " "; grep -F $i all_40peptide_results.tsv | cut -f10 | sort -u; done | sort -k2 > peptide_replacements.tsv
#manually fix
sed "s/  /\t/g" peptide_replacements.tsv | sort -k2 > test && mv test peptide_replacements.tsv

sed "s/ /_/g" all_40peptide_results.tsv | sort -k10 | join -1 10 -2 2 - peptide_replacements.tsv | sed "s/ /\t/g" | cut -f2- > all_40peptide_results.fixed.tsv
cut -f10,17 all_40peptide_results.fixed.tsv | sort -u | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/ //g" > peptide_protein.matches.tsv

cat all_40peptide_results.fixed.tsv | awk -F '\t' '($15<0.01)' | cut -f1,17 | rev | cut -f1 -d ":" | rev | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.01"}' > protein_distribution.csv
cat all_40peptide_results.fixed.tsv | awk -F '\t' '($15<0.001)' | cut -f1,17 | rev | cut -f1 -d ":" | rev | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.001"}' >> protein_distribution.csv
cat all_40peptide_results.fixed.tsv | awk -F '\t' '($15<0.0001)' | cut -f1,17 | rev | cut -f1 -d ":" | rev | cut -f2 | sort | uniq -c | awk '{print $2","$1",0.0001"}' >> protein_distribution.csv
sort -t ',' -k1 protein_distribution.csv -o protein_distribution.csv

sed "s/,/@/g" peptide_protein.matches.tsv | sed "s/\t/,/g" | sort -k1 -t ',' | join -t ',' -1 1 -2 1 - protein_distribution.csv | sed "1s/^/peptide,protein,conditions,criteria\n/g" | head > protein_PSM_distribution.csv

#R code:

setwd("/stor/scratch/Ochman/hassan/112724_protogene_extension/consolidating_datasets")

library(ggplot2)
library(dplyr)

# Read the data
p1 <- read.csv("protein_PSM_distribution.csv")

# Convert 'criteria' to factor for ordered layering
p1$criteria <- factor(p1$criteria, levels = c("0.0001q", "0.001q", "0.01q"))

# Split data for each criteria
p1_01 <- p1 %>% filter(criteria == "0.01q")
p1_001 <- p1 %>% filter(criteria == "0.001q")
p1_0001 <- p1 %>% filter(criteria == "0.0001q")

# Extract protein order from p1_01 based on decreasing "conditions" values
protein_order <- p1_01 %>% arrange(desc(conditions)) %>% pull(protein)

# Re-factor proteins in all datasets based on the extracted order
p1_01$protein <- factor(p1_01$protein, levels = protein_order)
p1_001$protein <- factor(p1_001$protein, levels = protein_order)
p1_0001$protein <- factor(p1_0001$protein, levels = protein_order)

# Create the overlaid bar plot
ggplot() +
  # Base layer: Criteria 0.01 (full bar height, lightest blue)
  geom_col(data = p1_01, aes(x = protein, y = conditions, fill = "0.01q"), width = 0.8) +
  # Middle layer: Criteria 0.001 (medium blue)
  geom_col(data = p1_001, aes(x = protein, y = conditions, fill = "0.001q"), width = 0.8) +
  # Top layer: Criteria 0.0001 (darkest blue)
  geom_col(data = p1_0001, aes(x = protein, y = conditions, fill = "0.0001q"), width = 0.8) +
  # Custom color scale
  scale_fill_manual(values = c("0.01q" = "#B3CDE3",  # Lightest blue
                               "0.001q" = "#6497B1", # Medium blue
                               "0.0001q" = "#005B96")) + # Darkest blue
  labs(x = "Protein", y = "Conditions", fill = "Criteria") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),  # Remove x-axis labels
        axis.ticks.x = element_line(color = "black"),
        axis.ticks.y = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major.x = element_blank(),  # Remove major vertical gridlines
        panel.grid.minor.x = element_blank(),
        legend.position = "none") +  # Remove legend
  scale_y_continuous(breaks = seq(0, 65, 5), limits = c(0, 65)) +  # Set y-axis intervals
  labs(x = "Proteins", y = "Number of conditions in which peptide is detected")  # Set axis labels

#manually fixed the case where one protein is supported by two distinct peptides

################

#Protein-peptide stuff

for i in $(grep -v "NZ_QOWW01000007.1_4686" peptide_protein.matches.tsv | cut -f2 | cut -f1 -d "," | sort -u)
do
grep -A1 "$i" all_qval0.01.prot.faa | seqkit fx2tab > temp
pept=$(grep "$i" peptide_protein.matches.tsv | cut -f1 | sed "s/+15.995//g" | sed "s/+57.021//g" | grep -o -F -f - temp)
echo "$i" > temp2
awk '{print length($2)}' temp >> temp2
grep "$i" peptide_protein.matches.tsv | cut -f1 | sed "s/+15.995//g" | sed "s/+57.021//g" | grep -F -f - temp | sed "s/$pept/@/g" | cut -f2 | sed "s/@/\t/g" | awk -F '\t' '{print length($1),length($2)}' >> temp2
sed "s/ /\n/g" temp2 | sed -z "s/\n/\t/g" | awk '{print $1"\t"1"\t"$3+1"\t"($2-$4)"\t"$2}' >> all_peptide_mapping.tsv
done
