#Expression:

#Map all Ecoli genes and protogenes to REL606:
cd /stor/work/Ochman/hassan/protogene_extension/Ecoli_list
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d REL606 -f 2 --gff3-fasta-annotation=1 Ecoli_CDS_queryfile.faa > Ecoli_CDS_queryfile.mapped.REL606.gff3
awk -F '\t' '($3=="mRNA")' Ecoli_CDS_queryfile.mapped.REL606.gff3 | grep "mrna1" | cut -f1 -d ";" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' | sed "s/.mrna1//g" | sed "s/ID=//g" > Ecoli_CDS_queryfile.mapped.REL606.gtf

awk -F '\t' '{OFS=FS}($1="REL606")' /stor/scratch/Ochman/hassan/112724_protogene_extension/REL606.getorf.bacterial.longest.final.gtf | bedtools intersect -v -s -a - -b Ecoli_CDS_queryfile.mapped.REL606.gtf | sed "s/NC_012967.1/control/g" > REL606.controlORFs.gtf
cat REL606.controlORFs.gtf Ecoli_CDS_queryfile.mapped.REL606.gtf > REL606_tpmcalculation.gtf

awk -F '\t' '{OFS=""}{print "time htseq-count -f bam -a 0 -t CDS --nonunique all /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/RNAseq\/bamfiles\/",$1,"_sorted_pairedonly_byname.bam REL606_tpmcalculation.gtf \| head -n -5 \| sed \"s\/^\/",$1,"\t\/g\" >> ",$1,"_htseq.tsv"}' /stor/work/Ochman/hassan/proteomics_denovo/1031_RNAseq_databasemaking/MURIfiles_with_RNAprotein_data.txt > running.sh

#Above step still running
#Make #windows above TPM per dataset across all conditions bar plot

awk -F '\t' '{print $5-$4,$9}' REL606_tpmcalculation.gtf | cut -f 1 -d ';' | sed "s/transcript id \"//g" | sed "s/\"//g" | sed "s/ transcript_id /\t/g" | awk -F '\t' '{OFS=FS}{print $2,$1}' | sort -k1 > ORF_lengths.txt

for i in MURI*htseq.tsv
do
sort -k2 $i | join -1 2 -2 1 - ORF_lengths.txt | sed "s/ /\t/g" | awk -F '\t' '{OFS=FS}{print $1, $2, $3, $4, $3/$4}' > interim
sum=$(awk '{ sum += $5 } END { print sum }' interim)
awk -v sum="$sum" 'BEGIN{OFS="\t"} {print $0, ($5/sum) * 1000000}' interim > $i.tpm
done

cat *tpm | cut -f1,2,6 > REL606_tpmcalculation.tpm

#Non-redundantify the REL606 annotated genes:

grep "REL606" ../Ecoli_list/Ecoli_annotated.gtf | awk '
{
    stop = ($7 == "+") ? $5 : $4  # Determine stop codon position
    feature_length = $5 - $4  # Compute feature length

    key = $1 ":" stop  # Unique key based on chromosome and stop position

    if (!(key in data) || feature_length > data[key]) {
        data[key] = feature_length
        lines[key] = $0  # Store the line corresponding to the longest entry
    }
}
END {
    for (k in lines) {
        print lines[k]
    }
}' > Ecoli_nonredundant_REL606.gtf


egrep -iv "prodigal|smorf|balrog|gms2" Ecoli_CDS_queryfile.mapped.REL606.gtf | egrep -v "NZ_" | bedtools intersect -s -wo -a - -b Ecoli_nonredundant_REL606.gtf | awk -F '\t' '($19>10)' | cut -f9 | sort -u | cut -f2 -d "\"" | sort -u > novelgenes_exclude.txt

egrep -iv "gms2|balrog|prodigal|smorf|NZ_" REL606_tpmcalculation.tpm | grep -v -F -f novelgenes_exclude.txt - | awk -F '_' '{print $1"\t"$0}' | sed "s/^EC/EC_Ndah/g" | sed "s/^NC/MS/g" | grep -v "control" | sed "s/$/\tnovel/g" > Ecoli.REL606.tpm
grep "control" REL606_tpmcalculation.tpm | sed "s/^/control\t/g" | sed "s/$/\tcontrol/g" >> Ecoli.REL606.tpm
egrep "gms2|balrog|prodigal|smorf" Ecoli_CDS_queryfile.mapped.REL606.tpm | sed "s/^/annotated\t/g" | sed "s/$/\tannotated/g" >> Ecoli.REL606.tpm
sed -i "1s/^/dataset\tgene\tcondition\ttpm\tcategory\n/g" Ecoli.REL606.tpm

awk -F '\t' '($4>0.3)' Ecoli.REL606.tpm | cut -f1,2,3,5 | tail -n+2 | sort -k3 | join -1 3 -2 1 - ../../proteomics_denovo/03052024_presenceabsence/MURI_uniqueID_tidy.tsv | awk '{print $2"@"$3"@"$4"\t"$5}' | sort -u  | cut -f1 | sort | uniq -c | awk '{print $2,$1}' | sed "s/ /\t/g" | sort -nrk2 | sed "s/@/\t/g" | sed "1s/^/dataset\tgene\tcategory\tconditions\n/g" > Ecoli.REL606.0.3tpm.conditions.tsv
awk -F '\t' '($4>1)' Ecoli.REL606.tpm | cut -f1,2,3,5 | tail -n+2 | sort -k3 | join -1 3 -2 1 - ../../proteomics_denovo/03052024_presenceabsence/MURI_uniqueID_tidy.tsv | awk '{print $2"@"$3"@"$4"\t"$5}' | sort -u  | cut -f1 | sort | uniq -c | awk '{print $2,$1}' | sed "s/ /\t/g" | sort -nrk2 | sed "s/@/\t/g" | sed "1s/^/dataset\tgene\tcategory\tconditions\n/g" > Ecoli.REL606.1tpm.conditions.tsv

cat Ecoli.REL606.tpm | cut -f1,2 | sort -u | tail -n+2 | cut -f1 | sort | uniq -c | head -n-1 | awk '{print $2"\t"$1}' | sed "s/$/\ttotal/g" > Ecoli.REL606.tpm.numbers
awk -F '\t' '($4>1)' Ecoli.REL606.tpm | cut -f1,2 | sort -u | tail -n+2 | cut -f1 | sort | uniq -c | head -n-1| awk '{print $2"\t"$1}' | sed "s/$/\t1tpm/g" >> Ecoli.REL606.tpm.numbers
awk -F '\t' '($4>5)' Ecoli.REL606.tpm | cut -f1,2 | sort -u | tail -n+2 | cut -f1 | sort | uniq -c | head -n-1| awk '{print $2"\t"$1}' | sed "s/$/\t5tpm/g" >> Ecoli.REL606.tpm.numbers
sed -i "1s/^/dataset\tcount\tcategory\n/g" Ecoli.REL606.tpm.numbers

#R codes
#Cumdist:

p1 <- read.csv("Ecoli.REL606.1tpm.conditions.tsv",sep='\t')
p2 <- read.csv("Ecoli.REL606.0.3tpm.conditions.tsv",sep='\t')

generate_plot_data <- function(df) {
  cutoffs <- seq(-0.001, 41, by = 1)
  plot_data_list <- list()
  
  for (line_value in unique(df$dataset)) {
    line_data <- df[df$dataset == line_value, ]
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
line_colors <- c(
  "EC_Ndah" = "darkmagenta",
  "MS" = "magenta",
  "Nakahigashi" = "pink",
  "Stringer" = "blue",
  "VanOrsdel" = "green",
  "Weaver" = "azure4",
  "annotated" = "orange"
)

# Create ggplot cumulative distribution plot for p1
plot_p1 <- ggplot(plot_data_p1, aes(x = Cutoff, y = Fraction_Passed, color = Line)) +
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
plot_p2 <- ggplot(plot_data_p2, aes(x = Cutoff, y = Fraction_Passed, color = Line)) +
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
combined_plot <- plot_p1 + plot_p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Print the combined plot
print(combined_plot)

###TPM count###

library(ggplot2)
library(dplyr)

p1 <- read.csv("Ecoli.REL606.tpm.numbers",sep='\t')

dataset_order <- c("annotated", "Nakahigashi", "Stringer", "Weaver", 
                   "VanOrsdel", "EC_Ndah", "MS", "control")

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
