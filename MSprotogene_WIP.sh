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

