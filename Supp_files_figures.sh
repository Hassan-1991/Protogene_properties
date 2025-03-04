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

ls /stor/scratch/Ochman/hassan/112724_protogene_extension/promising_spectra/*spectra.tsv
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

######Expression analyses, as requested by reviewr 1######

#Questions:

1. Does expression vary depending on conservation? For either annotated or novel proteins
    a) Cumdist
    b) Mean TPM (violin plot) of groups across all conditions
    c) What fraction passes TPM cutoffs (1/5 TPM)
2. Is there any correlation between being consistently translated ((i)detected across multiple datasets, or (ii) MS detected across lots of conditions) and level/consistency of expression
    a) Cumdist
    b) Mean TPM (violin plot) of groups across all conditions
    c) What fraction passes TPM cutoffs (1/5 TPM)
3. Can any difference be explained in terms of readthrough?

#1(a) - Cumdist

grep "control" Ecoli.REL606.1tpm.conditions.tsv | awk '{OFS=FS}{print $2"\t"$1"\t"$4}' >> Ecoli.REL606.1tpm.conditions.conservation.tsv
grep "control" Ecoli.REL606.0.3tpm.conditions.tsv | awk '{OFS=FS}{print $2"\t"$1"\t"$4}' >> Ecoli.REL606.0.3tpm.conditions.conservation.tsv

p1 <- read.csv("Ecoli.REL606.1tpm.conditions.conservation.tsv",sep='\t')
p2 <- read.csv("Ecoli.REL606.0.3tpm.conditions.conservation.tsv",sep='\t')

generate_plot_data <- function(df) {
  cutoffs <- seq(-0.001, 41, by = 1)
  plot_data_list <- list()
  
  for (line_value in unique(df$category)) {
    line_data <- df[df$category == line_value, ]
    fraction_passed <- numeric(length(cutoffs))
    
    for (i in seq_along(cutoffs)) {
      fraction_passed[i] <- sum(line_data$conditions > cutoffs[i]) / nrow(line_data)
    }
    
    plot_data <- data.frame(Cutoff = cutoffs, Fraction_Passed = fraction_passed, Line = line_value)
    plot_data_list[[line_value]] <- plot_data
  }
  
  combined_plot_data <- do.call(rbind, plot_data_list)
  
  return(combined_plot_data)
}

# Generate plot data for p1 and p2
plot_data_p1 <- generate_plot_data(p1)
plot_data_p2 <- generate_plot_data(p2)

# Define line colors with new labels
#line_colors <- c(
#  "EC_Ndah" = "",
#  "MS" = "magenta",
#  "Nakahigashi" = "#4f83cc",
#  "Stringer" = "#8da0cb",
#  "VanOrsdel" = "green",
#  "Weaver" = "#08306b",
#  "annotated" = "darkmagenta",
#  "control" = "grey"
#)

line_colors <- c(
  "annotated_conserved" = "darkmagenta",
  "annotated_ORFan" = "magenta",
  "novel_conserved" = "#08306b",
  "novel_ORFan" = "#8da0cb",
  "control" = "grey"
)

# Create ggplot cumulative distribution plot for p1
plot_p1 <- ggplot(plot_data_p1 %>% filter(Line!="MS") %>% filter(Line!="EC_Ndah") %>% filter(Line!="VanOrsdel"), aes(x = Cutoff, y = Fraction_Passed, color = Line)) +
  geom_line(size = 1) +
  labs(
    x = "Number of unique conditions in which genes are expressed",
    y = "Fraction of genes in each category"
  ) +
  scale_x_continuous(breaks = seq(0, 42, by = 3)) +
  scale_color_manual(values = line_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20),  # Adjust the size as needed
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24), # Adjust the size as needed
    axis.title.y = element_text(size = 24), # Adjust the size as needed
    legend.text = element_text(size = 12),  # Adjust the size as needed
    legend.title = element_text(size = 14),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black", size = 0.7)
  )

# Create ggplot cumulative distribution plot for p2
plot_p2 <- ggplot(plot_data_p2 %>% filter(Line!="MS") %>% filter(Line!="EC_Ndah") %>% filter(Line!="VanOrsdel"), aes(x = Cutoff, y = Fraction_Passed, color = Line)) +
  geom_line(size = 1) +
  labs(
    x = "Number of unique conditions in which genes are expressed",
    y = "Fraction of genes in each category"
  ) +
  scale_x_continuous(breaks = seq(0, 42, by = 3)) +
  scale_color_manual(values = line_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20),  # Adjust the size as needed
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 24), # Adjust the size as needed
    axis.title.y = element_text(size = 24), # Adjust the size as needed
    legend.position = "none",  # Remove the legend from the second plot
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.ticks = element_line(color = "black", size = 0.7)
  )

# Combine the two plots side-by-side with a shared legend
cumdist_combined_plot <- plot_p1 + plot_p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Print the combined plot
print(cumdist_combined_plot)

#1(b) - MeanTPM

awk -F'\t' '{sum[$1,$4] += $3; count[$1,$4]++; category[$1] = $4} 
     END {for (key in sum) {split(key, arr, SUBSEP); print arr[1], arr[2], sum[key] / count[key]}}' OFS='\t' Ecoli.REL606.tpm.conservation.tsv | sed "1s/^/gene\tcategory\tmeantpm\n/g" > Ecoli.REL606.tpm.conservation.meanvalues.tsv

#1(c) - TPM cutoff

setwd("/stor/work/Ochman/hassan/protogene_extension/expression_location_properties")

p1 <- read.csv("Ecoli.REL606.tpm.numbers.conservation",sep='\t')

dataset_order <- c("annotated_conserved", "annotated_ORFan", "novel_conserved", "novel_ORFan")

# Normalize counts relative to the 'total' values
p1_normalized <- p1 %>%
  group_by(dataset) %>%
  mutate(fraction = ifelse(category != "total", count / count[category == "total"], NA)) %>%
  filter(category %in% c("1tpm", "5tpm"))  # Exclude the "total" rows from plotting

p1_normalized$dataset <- factor(p1_normalized$dataset, levels = dataset_order)

# Create grouped bar plot
ggplot(p1_normalized, aes(x = dataset, y = fraction, fill = category)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(x = "Dataset", y = "Fraction", fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


ggplot(p1_normalized, aes(x = dataset, y = fraction, fill = category)) +
  geom_col(data = p1_normalized %>% filter(category == "1tpm"), position = "identity", fill = "steelblue") +  # Full 1tpm bars
  geom_col(data = p1_normalized %>% filter(category == "5tpm"), position = "identity", fill = "orange") +  # Smaller 5tpm bars
  theme_minimal() +
  labs(x = "Dataset", y = "Fraction", fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggplot(p1_normalized, aes(x = dataset, y = fraction, fill = category)) +
  geom_col(data = p1_normalized %>% filter(category == "1tpm"), position = "identity", 
           fill = "lightblue", color = "black", width = 0.7) +  # Light blue with black border
  geom_col(data = p1_normalized %>% filter(category == "5tpm"), position = "identity", 
           fill = "dodgerblue4", color = "black", width = 0.7) +  # Darker blue with black border
  theme_minimal(base_size = 14) +
  labs(x = "Dataset", y = "Fraction", fill = "Category") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),  # Box around the plot
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    axis.ticks = element_line(color = "black"),  # Add tick marks
    axis.line.x = element_line(color = "black", size = 0.25),  # Line along the base of bars
    panel.grid.major.y = element_line(color = "gray80", linetype = "solid"),  # Horizontal gridlines
    panel.grid.minor = element_blank()  # Keep minor gridlines off
  )

#Now consistency analysis, including expression

cd Ecoli_list
grep "^>" Salmonella_protein_queryfile.faa | egrep -iv "prodigal|smorf|balrog|gms2" | grep -vf ../expression_location_properties/sequence_properties/protogene_exclude.txt - | grep --no-group-separator -A1 -f - Salmonella_CDS_queryfile.faa > consistency_analysis/Salmonella_135_protogenes.CDS.faa
/stor/work/Ochman/hassan/tools/faTrans -stop Salmonella_135_protogenes.CDS.faa Salmonella_135_protogenes.prot.faa
seqkit fx2tab Salmonella_135_protogenes.prot.faa | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > temp && mv temp Salmonella_135_protogenes.prot.faa
#Same w Ecoli

usearch -sortbylength Ecoli_559_protogenes.prot.faa -fastaout Ecoli_559_protogenes.prot.sorted.faa -minseqlength 1
usearch -cluster_smallmem Ecoli_559_protogenes.prot.sorted.faa -id 0.9 -centroids Ecoli_559_protogenes.prot.sorted.nr.faa -uc Ecoli_559_protogenes.prot.clusters.uc
awk '$1=="S"{c[$2]=$9}$1=="H"{h[$2]=(h[$2]?h[$2]","$9:$9)}END{for(i in h)print "Cluster " i ":\t" c[i] "," h[i]}' Ecoli_559_protogenes.prot.clusters.uc > Ecoli_consistency.interim

usearch -sortbylength Salmonella_135_protogenes.prot.faa -fastaout Salmonella_135_protogenes.prot.sorted.faa -minseqlength 1
usearch -cluster_smallmem Salmonella_135_protogenes.prot.sorted.faa -id 0.9 -centroids Salmonella_135_protogenes.prot.sorted.nr.faa -uc Salmonella_135_protogenes.prot.clusters.uc
awk '$1=="S"{c[$2]=$9}$1=="H"{h[$2]=(h[$2]?h[$2]","$9:$9)}END{for(i in h)print "Cluster " i ":\t" c[i] "," h[i]}' Salmonella_135_protogenes.prot.clusters.uc > Salmonella_consistency.interim

#Manually parse to remove hits within the same database

#22 are consistent between databases - pretty high number!
#Ndah and Fijalkowski probably overlaps...

#37 are consistent
#How many are robustly translated?

cut -f1 -d ',' Ecoli_consistency.interim | cut -f2 | sort -u | grep -F -f - Ecoli_559_protogenes.prot.sorted.nr.faa | tr -d ">" | grep -F -f - ../../expression_location_properties/Ecoli.REL606.1tpm.conditions.tsv | wc -l
#i.e., 28 are independent of annotated transcripts
#Of these 28, 25 are expressed in >10 unique conditions:
cut -f1 -d ',' Ecoli_consistency.interim | cut -f2 | sort -u | grep -F -f - Ecoli_559_protogenes.prot.sorted.nr.faa | tr -d ">" | grep -F -f - ../../expression_location_properties/Ecoli.REL606.1tpm.conditions.tsv | awk -F '\t' '($4>10)' | wc -l
#As opposed to 195 of all proto-genes:
egrep -v "control|smorf|prod|gms2|balrog" ../../expression_location_properties/Ecoli.REL606.1tpm.conditions.tsv | awk -F '\t' '($4>10)' | wc -l
#Which is (170/309=) 63% of total:
egrep -v "control|smorf|prod|gms2|balrog" ../../expression_location_properties/Ecoli.REL606.1tpm.conditions.tsv | wc -l

#How many are ORFans? 18:
cut -f1 -d ',' Ecoli_consistency.interim | cut -f2 | sort -u | grep -F -f - Ecoli_559_protogenes.prot.sorted.nr.faa | tr -d ">" | grep -F -f - ../../expression_location_properties/sequence_properties/Ecoli_genusspecific_ORFans.final.txt | wc -l
#I.e., 18/37 =  about half of these are ORFans, which is comparable to the rate of novel genes in general. Nothing interesting follows from this/.l

#Conservation

egrep "prodigal|balrog|smorf|gms2" Mycobacterium_extragenus_hits.genusnumber | awk '{print $1}' | sort -n | awk 'NR%2==1{a[NR]=$0} END{if(NR%2==1)print a[int(NR/2)+1]; else print (a[NR/2]+a[NR/2+1])/2}'
#181.5
egrep -iv "prodigal|balrog|smorf|gms2" Mycobacterium_extragenus_hits.genusnumber | awk '{print $1}' | sort -n | awk 'NR%2==1{a[NR]=$0} END{if(NR%2==1)print a[int(NR/2)+1]; else print (a[NR/2]+a[NR/2+1])/2}'
#4

#For Ecoli, the values are 159.5 and 6
#For Salmnoella, however, the values are 618 and 18

#LET'S MAKE THE LARGE TABLE#

mkdir /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/supp_allprotein_categories


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
