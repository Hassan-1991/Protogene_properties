#This code curates all repported proto-genes from E. coli, Salmonella enterica and M. tuberculosis

#Studies chosen for E. coli:

#Nakahigashi_2016
#Ndah_2017
#VanOrsdel_2018
#Weaver_2019
#Stringer_2021

#Studies chosen for S. enterica:

#Baek_2017
#Giess_2017
#Ndah_2017
#Venturini_2020
#Willems_2020
#Fijalkowski_2022

#Study chosen for M. tuberculosis:

#Smith_2022

#These studies provide the genome coordinates of novel proteins
#The task is to identify what genomes these coordinates are from and extract
#After some trial-and-error, here are the genomes from which the sequences are to be extracted:

#Stringer_2021 - CP001509_3.fna (strain BL21)
#Weaver_2019, VanOrsdel_2018, Ndah_2017 - sequence_RS.fasta (updated version of MG1655)
#Nakahigashi_2016 - sequence_oldMG1655.fasta (old version of MG1655)
#My MS-detected cases - REL606.fasta

#Baek_2017 - Salmonella_14028s_Baek_2017.fasta
#Giess_2017, Ndah_2017, Venturini_2020, Willems_2020, Fijalkowski_2022 - GCA_000210855.2_ASM21085v2_genomic.fna

#Coordinates are variations of supplemental files, so different approaches are used to extract the info

#Ecoli: Stringer_2021: Only retain sequences greater or equal to 30bp
awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$3-$2}' Stringer_2021 | tail -n+2 | grep -v "^y" | sed "s/-\t-/-\t/g" | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+1}' | awk -F '\t' '($5>=30)' > Stringer_2021_filtered
awk -F '\t' '{print "CP001509.3\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Stringer_"NR"\";gene_id \"Stringer_"NR"\";"}' Stringer_2021_filtered | awk -F '\t' '($7=="+")' > Stringer_2021.gtf
awk -F '\t' '{print "CP001509.3\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Stringer_"NR"\";gene_id \"Stringer_"NR"\";"}' Stringer_2021_filtered | awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$4,$6,$7,$8,$9}' >> Stringer_2021.gtf

#Ecoli: Weaver_2019: Only retain those with amino acid lengths more than 9
awk -F '\t' '{OFS=FS}{print $1,$2,$9+1,$10}' Weaver_2019 | tail -n+2 | sed "s/plus/+/g" | sed "s/minus/-/g" | awk -F '\t' '($3>9)' > Weaver_2019_filtered
awk -F '\t' '{print "NC_000913.3\t.\tCDS\t"$1"\t"$2"\t.\t"$4"\t.\ttranscript_id \"Weaver_"NR"\";gene_id \"Weaver_"NR"\";"}' Weaver_2019_filtered > Weaver_2019.gtf

#Ecoli: VanOrsdel_2018: Only retain those that get a "yes" on the expression column
awk -F '\t' '($2=="Y")' VanOrsdel_2018 | cut -f 1,6 | sed "s/, /\t/g" > VanOrsdel_2018_filtered
awk -F '\t' '{print "NC_000913.3\t.\tCDS\t"$2"\t"$3"\t.\t+\t.\ttranscript_id \""$1"\";gene_id \""$1"\";"}' VanOrsdel_2018_filtered | sed "s/ //g" | sed "s/_id/_id /g" | gtf2bed | bedtools getfasta -s -name -fi sequence_RS.fasta -bed - | sed '/^>/ s/:.*//' | seqkit fx2tab | egrep $'\tATG|\tGTG' | rev | sed "s/\t//" | rev | sed "s/$/#/g" | egrep "TGA#|TAA#" | cut -f1 -d "#" | cut -f1 | grep -f - VanOrsdel_2018_filtered | sed "s/$/\t+/g" > VanOrsdel_2018_filtered_final
awk -F '\t' '{print "NC_000913.3\t.\tCDS\t"$2"\t"$3"\t.\t-\t.\ttranscript_id \""$1"\";gene_id \""$1"\";"}' VanOrsdel_2018_filtered | sed "s/ //g" | sed "s/_id/_id /g" | gtf2bed | bedtools getfasta -s -name -fi sequence_RS.fasta -bed - | sed '/^>/ s/:.*//' | seqkit fx2tab | egrep $'\tATG|\tGTG' | rev | sed "s/\t//" | rev | sed "s/$/#/g" | egrep "TGA#|TAA#" | sed "s/id i/idi/g" | cut -f1 | grep -f - VanOrsdel_2018_filtered | sed "s/$/\t-/g" >> VanOrsdel_2018_filtered_final
awk -F '\t' '{print "NC_000913.3\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"VanOrsdel_"NR"\";gene_id \"VanOrsdel_"NR"\";"}' VanOrsdel_2018_filtered_final > VanOrsdel_2018.gtf

#Ecoli: Ndah_2017
tail -n+2 Ndah_2017 | cut -f1-3 | cut -f 2- -d ":" | sed "s/\t-\t/\t%\t/g" | sed "s/-/\t/g" | sed "s/%/-/g" > Ndah_2017_filtered
awk -F '\t' '{print "NC_000913.3\t.\tCDS\t"$1"\t"$2"\t.\t"$3"\t.\ttranscript_id \"Ndah_"NR"\";gene_id \"Ndah_"NR"\";"}' Ndah_2017_filtered | awk -F '\t' '($7=="+")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+3,$6,$7,$8,$9}' > Ndah_2017.gtf
awk -F '\t' '{print "NC_000913.3\t.\tCDS\t"$1"\t"$2"\t.\t"$3"\t.\ttranscript_id \"Ndah_"NR"\";gene_id \"Ndah_"NR"\";"}' Ndah_2017_filtered | awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-3,$5,$6,$7,$8,$9}' >> Ndah_2017.gtf

#Ecoli: Nakahigashi_2016: Only retain those that get a "yes" on the expression column
awk -F '\t' '($9=="sORF")' Nakahigashi_2016 | awk -F '\t' '{OFS=FS}{print $4,$6,$7,$8+1}' | awk -F '\t' '($4>9)' > Nakahigashi_2016_filtered
awk -F '\t' '{print "U00096.2\t.\tCDS\t"$2"\t"$3"\t.\t"$1"\t.\ttranscript_id \"Nakahigashi_"NR"\";gene_id \"Nakahigashi_"NR"\";"}' Nakahigashi_2016_filtered | awk -F '\t' '($7=="+")' > Nakahigashi_2016.gtf
awk -F '\t' '{print "U00096.2\t.\tCDS\t"$2"\t"$3"\t.\t"$1"\t.\ttranscript_id \"Nakahigashi_"NR"\";gene_id \"Nakahigashi_"NR"\";"}' Nakahigashi_2016_filtered | awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$5,$4,$6,$7,$8,$9}' >> Nakahigashi_2016.gtf

#Ecoli: My mass spec proteins from five strains:
#REL606:
grep "^>" /stor/scratch/Ochman/hassan/112724_protogene_extension/Caglar2017_MSvalidated_proteins.cds.faa | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - /stor/scratch/Ochman/hassan/112724_protogene_extension/REL606.final.gtf | awk -F '\t' '{OFS=FS}{print "REL606",$2,$3,$4,$5,$6,$7,$8,$9}' > MS_REL606.gtf
#K12MG1655:
sed "s/.*/\"&\"/g" /stor/scratch/Ochman/hassan/112724_protogene_extension/Mori2021 | grep -F -f - /stor/scratch/Ochman/hassan/112724_protogene_extension/K12MG1655.final.gtf > MS_K12MG1655.gtf
#ECOR:
sed "s/.*/\"&\"/g" /stor/scratch/Ochman/hassan/112724_protogene_extension/ECOR2023_MSvalidated_proteins.txt | grep -F -f - /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/*genome.final.gtf | grep "ECOR_11" | cut -f2- -d ":" > MS_ECOR_11.gtf
sed "s/.*/\"&\"/g" /stor/scratch/Ochman/hassan/112724_protogene_extension/ECOR2023_MSvalidated_proteins.txt | grep -F -f - /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/*genome.final.gtf | grep "ECOR_27" | cut -f2- -d ":" > MS_ECOR_27.gtf
sed "s/.*/\"&\"/g" /stor/scratch/Ochman/hassan/112724_protogene_extension/ECOR2023_MSvalidated_proteins.txt | grep -F -f - /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/*genome.final.gtf | grep "ECOR_37" | cut -f2- -d ":" > MS_ECOR_37.gtf

#Salmonella: Fijalkowski_2022

awk -F '\t' '($5=="Intergenic")' Fijalkowski_2022 | cut -f2,3 | sed "s/Chromosome/FQ312003.1/g" | sed "s/pCol1B9_SL1344/HE654725.1/g" | sed "s/pSLT_SL1344/HE654724.1/g" | sed "s/:/\t/g" | sed "s/\t-/\t%/g" | sed "s/-/\t/g" | sed "s/%/-/g" | awk -F '\t' '{print $1"\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Fijalkowski_"NR"\";gene_id \"Fijalkowski_"NR"\";"}' | sed "s/52360\t52034/52034\t52360/g" | awk -F '\t' '($7=="+")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+3,$6,$7,$8,$9}' > Fijalkowski_2022.gtf
awk -F '\t' '($5=="Intergenic")' Fijalkowski_2022 | cut -f2,3 | sed "s/Chromosome/FQ312003.1/g" | sed "s/pCol1B9_SL1344/HE654725.1/g" | sed "s/pSLT_SL1344/HE654724.1/g" | sed "s/:/\t/g" | sed "s/\t-/\t%/g" | sed "s/-/\t/g" | sed "s/%/-/g" | awk -F '\t' '{print $1"\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Fijalkowski_"NR"\";gene_id \"Fijalkowski_"NR"\";"}' | sed "s/52360\t52034/52034\t52360/g" | awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-3,$5,$6,$7,$8,$9}' >> Fijalkowski_2022.gtf

#Salmonella: Willems_2020
awk -F '\t' '($3=="Intergenic"||$3=="ncRNA"||$3=="Indel"||$3=="Frameshift")' Willems_2020 | cut -f1 | sed "s/Chromosome/FQ312003.1/g" | sed "s/pCol1B9_SL1344/HE654725.1/g" | sed "s/pSLT_SL1344/HE654724.1/g" | sed "s/:/\t/g" | sed "s/-/\t/" | sed "s/_/\t/g" | awk -F '\t' '{print $1"\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Willems_"NR"\";gene_id \"Willems_"NR"\";"}' | awk -F '\t' '($7=="+")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+3,$6,$7,$8,$9}' > Willems_2020.gtf
awk -F '\t' '($3=="Intergenic"||$3=="ncRNA"||$3=="Indel"||$3=="Frameshift")' Willems_2020 | cut -f1 | sed "s/Chromosome/FQ312003.1/g" | sed "s/pCol1B9_SL1344/HE654725.1/g" | sed "s/pSLT_SL1344/HE654724.1/g" | sed "s/:/\t/g" | sed "s/-/\t/" | sed "s/_/\t/g" | awk -F '\t' '{print $1"\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Willems_"NR"\";gene_id \"Willems_"NR"\";"}' | awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-3,$5,$6,$7,$8,$9}' >> Willems_2020.gtf

#Salmonella: Venturini_2020
awk -F '\t' '($6=="y")' Venturini_2020 | cut -f 2-4 | awk -F '\t' '{print "FQ312003.1\t.\tCDS\t"$1"\t"$2"\t.\t"$3"\t.\ttranscript_id \"Venturini_"NR"\";gene_id \"Venturini_"NR"\";"}' > Venturini_2020.gtf

#Salmonella: Ndah_2017
tail -n+2 ST_Ndah_2017 | cut -f1,2 | sed "s/Chromosome/FQ312003.1/g" | sed "s/pCol1B9_SL1344/HE654725.1/g" | sed "s/pSLT_SL1344/HE654724.1/g" | sed "s/:/\t/g" | sed "s/-/\t/" | awk -F '\t' '{print $1"\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Ndah_"NR"\";gene_id \"Ndah_"NR"\";"}' | awk -F '\t' '($7=="+")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4,$5+3,$6,$7,$8,$9}' > ST_Ndah_2017.gtf
tail -n+2 ST_Ndah_2017 | cut -f1,2 | sed "s/Chromosome/FQ312003.1/g" | sed "s/pCol1B9_SL1344/HE654725.1/g" | sed "s/pSLT_SL1344/HE654724.1/g" | sed "s/:/\t/g" | sed "s/-/\t/" | awk -F '\t' '{print $1"\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Ndah_"NR"\";gene_id \"Ndah_"NR"\";"}' | awk -F '\t' '($7=="-")' | awk -F '\t' '{OFS=FS}{print $1,$2,$3,$4-3,$5,$6,$7,$8,$9}' >> ST_Ndah_2017.gtf

#Salmonella: Giess_2017
tail -n+2 Giess_2017 | cut -f1,2,3,4 | sed "s/Chromosome/FQ312003.1/g" | sed "s/pCol1B9_SL1344/HE654725.1/g" | sed "s/pSLT_SL1344/HE654724.1/g" | sed "s/pRSF1010_SL1344/HE654726.1/g" | awk -F '\t' '{print $1"\t.\tCDS\t"$2"\t"$3"\t.\t"$4"\t.\ttranscript_id \"Giess_"NR"\";gene_id \"Giess_"NR"\";"}' > Giess_2017.gtf

#Salmonella: Baek_2017
tail -n+2 Baek_2017  | cut -f 1,3,4 | grep -v "^A" | awk -F '\t' '{print "CP001363.1\t.\tCDS\t"$2"\t"$3"\t.\t+\t.\ttranscript_id \""$1"\";gene_id \""$1"\";"}' | sed "s/ (from pseudogene)//g" | gtf2bed | bedtools getfasta -s -name -fi Salmonella_14028s.fasta -bed - | seqkit fx2tab | sed "s/$/#/g" | rev | sed "s/\t//" | rev | sed "s/\t/%/g" | egrep "%ATG|%TTG|%GTG" | egrep "TAA#|TGA#|TAG#" | cut -f1 -d "(" | grep -w -f - Baek_2017 | cut -f 3,4 | grep -v "^M" | awk -F '\t' '{print "CP001363.1\t.\tCDS\t"$1"\t"$2"\t.\t+\t.\ttranscript_id \"Baekplus_"NR"\";gene_id \"Baekplus_"NR"\";"}' > Baek_2017.gtf
tail -n+2 Baek_2017  | cut -f 1,3,4 | grep -v "^A" | awk -F '\t' '{print "CP001363.1\t.\tCDS\t"$2"\t"$3"\t.\t-\t.\ttranscript_id \""$1"\";gene_id \""$1"\";"}' | gtf2bed | bedtools getfasta -s -name -fi Salmonella_14028s.fasta -bed - | seqkit fx2tab | sed "s/$/#/g" | rev | sed "s/\t//" | rev | sed "s/\t/%/g" | egrep "%ATG|%TTG|%GTG" | egrep "TAA#|TGA#|TAG#" | cut -f1 -d "(" | grep -w -f - Baek_2017 | cut -f 3,4 | grep -v "^M" | awk -F '\t' '{print "CP001363.1\t.\tCDS\t"$1"\t"$2"\t.\t-\t.\ttranscript_id \"Baekminus_"NR"\";gene_id \"Baekminus_"NR"\";"}' >> Baek_2017.gtf

#Put them all on a gtf
cat Stringer_2021.gtf Weaver_2019.gtf VanOrsdel_2018.gtf Ndah_2017.gtf Nakahigashi_2016.gtf MS_REL606.gtf MS_K12MG1655.gtf MS_ECOR_11.gtf MS_ECOR_27.gtf MS_ECOR_37.gtf | sed "s/Ndah/EC_Ndah/g" > all_protogenes.gtf
cat ../Salmonella_list/Fijalkowski_2022.gtf ../Salmonella_list/Willems_2020.gtf ../Salmonella_list/Venturini_2020.gtf ../Salmonella_list/Ndah_2017.gtf ../Salmonella_list/Giess_2017.gtf ../Salmonella_list/Baek_2017.gtf | sed "s/Ndah/ST_Ndah/g" >> all_protogenes.gtf
cat ../Mycobacterium_list/Mycobacterium_all_protogenes_final.gtf >> all_protogenes.gtf

#How to identify species-subsets of proto-genes:

#Mycobac:
egrep -i "leaderless|riboret" all_protogenes_tobeexcluded.txt
#Salmonella:
egrep -i "Giess|ST_Ndah|Venturini|Fijalkowski" all_protogenes_tobeexcluded.txt
#Ecoli:
egrep -i "Stringer|VanOrsdel|Weaver|Nakahigashi|EC_Ndah|NC_|NZ_" all_protogenes_tobeexcluded.txt

#Annotate ALL genomes using prodigal, genemarks2 and balrog

mkdir all_genomes
#List of genomes:
#Ecoli:
1. CP001509_3.fna
2. REL606
3. U00096.2
4. NC_000913.3
5, 6, 7. ECOR11,27,37 (balrog continuing)
#Balrog didn't work for one-quarter of ECOR11 and ECOR27
cat ECOR11_aa_balrog.gff ECOR11_ab_ab_balrog.gff > ECOR_11_genome_balrog.gff #incomplete
cat ECOR27_aa_ab_balrog.gff ECOR27_ab_balrog.gff > ECOR_27_genome_balrog.gff #incomplete
cat ECOR37_aa_balrog.gff ECOR37_ab_balrog.gff > ECOR_37_genome_balrog.gff

#Salmonella:
../Salmonella_list/GCA_000210855.2_ASM21085v2_genomic.fna
#CP001363.1
../Salmonella_list/Salmonella_14028s.fasta
#FQ312003.1
#HE654725.1
#HE654726.1
#HE654724.1

#Mycobacterium:
../Mycobacterium_list/MTb_novel_coordinates.gtf

#Finish annotation and curate the rest of these

#REL606, K12MG1655, Salmonella, 
cp /stor/work/Ochman/hassan/tools/gms2_linux_64/.gmhmmp2_key .
for i in sequence_RS.fasta sequence_oldMG1655.fasta CP001509_3.fna REL606.fasta Salmonella_14028s.fasta GCA_000210855.2_ASM21085v2_genomic.fna; do echo "/stor/work/Ochman/hassan/tools/prod_gms_balrog.sh ${i}"; done > running.sh
/stor/work/Ochman/hassan/tools/parallelize_run.sh

#Genemarks ain't running, troubleshoot
#annotate for small proteins:
conda activate smorfinder_env
for i in sequence_RS.fasta sequence_oldMG1655.fasta CP001509_3.fna REL606.fasta Salmonella_14028s.fasta GCA_000210855.2_ASM21085v2_genomic.fna; do smorf single $i && mv smorf_output "$i"_smorf_output; done

#Turn the gffs into gtfs

for i in ECOR_11_genome ECOR_27_genome ECOR_37_genome
do
cut -f-8 "$i"*balrog.gff | awk -F '\t' '{print $0,"\ttranscript_id \""$1"_balrog_"NR"\";gene_id \""$1"_balrog_"NR"\";"}' > "$i"_balrog.gtf
done

for i in sequence_RS sequence_oldMG1655 CP001509_3 REL606 Salmonella_14028s GCA_000210855.2
do
cut -f-8 "$i"*smorf_output/*gff | awk -F '\t' '{print $0,"transcript_id \""$1"_smorfer_"NR"\";gene_id \""$1"_smorfer_"NR"\";"}' > "$i"_smorfer.gtf
cut -f-8 "$i"*balrog.gff | awk -F '\t' '{print $0,"transcript_id \""$1"_balrog_"NR"\";gene_id \""$1"_balrog_"NR"\";"}' > "$i"_balrog.gtf
cut -f-8 "$i"*prodigal.gff | awk -F '\t' '{print $0,"transcript_id \""$1"_prodigal_"NR"\";gene_id \""$1"_prodigal_"NR"\";"}' > "$i"_prodigal.gtf
awk -F '\t' '($3=="CDS")' "$i"*genemarks2.gff | cut -f -8 | awk -F '\t' '{print $0,"transcript_id \""$1"_gms2_"NR"\";gene_id \""$1"_gms2_"NR"\";"}' > "$i"_gms2.gtf
done

#cat em
cat *smorfer.gtf *balrog.gtf *prodigal.gtf *gms2.gtf > annotated.gtf
cat ../Salmonella_list/*smorfer.gtf ../Salmonella_list/*balrog.gtf ../Salmonella_list/*prodigal.gtf ../Salmonella_list/*gms2.gtf >> annotated.gtf
cat ../Mycobacterium_list/*smorfer.gtf ../Mycobacterium_list/*balrog.gtf ../Mycobacterium_list/*prodigal.gtf ../Mycobacterium_list/*gms2.gtf >> annotated.gtf

bedtools intersect -wo -s -a all_protogenes.gtf -b annotated.gtf | awk -F '\t' '($7==$16&&$7=="+")' | awk -F '\t' '($5==$14)' | cut -f9 | cut -f 2 -d "\"" | sort -u > all_protogenes_tobeexcluded.txt
bedtools intersect -wo -s -a all_protogenes.gtf -b annotated.gtf | awk -F '\t' '($7==$16&&$7=="-")' | awk -F '\t' '($4==$13)' | cut -f9 | cut -f 2 -d "\"" | sort -u >> all_protogenes_tobeexcluded.txt

#Final proto-genes:

sed "s/.*/\"&\"/g" all_protogenes_tobeexcluded.txt | grep -vf- all_protogenes.gtf | egrep -i "leaderless|riboret" > Mycobacterium_protogenes.final.gtf
sed "s/.*/\"&\"/g" all_protogenes_tobeexcluded.txt | grep -vf- all_protogenes.gtf | egrep -i "Giess|ST_Ndah|Venturini|Fijalkowski" > Salmonella_protogenes.final.gtf
sed "s/.*/\"&\"/g" all_protogenes_tobeexcluded.txt | grep -vf- all_protogenes.gtf | egrep -i "Stringer|VanOrsdel|Weaver|Nakahigashi|EC_Ndah|\"NC_|\"NZ_" > Ecoli_protogenes.final.gtf

#Final annotated genes are non-redundant versions of genes annotated from the same strains whence the proto-genes came

egrep "^NZ_|CP001509.3|NC_000913.3|REL606|U00096.2" annotated.gtf > Ecoli_annotated.gtf
egrep "CP001363.1|FQ312003.1|HE654724.1|HE654725.1|HE654726.1" annotated.gtf > Salmonella_annotated.gtf
egrep "NC_99" annotated.gtf > Mycobacterium_annotated.gtf

cat CP001509_3.fna ECOR_11_genome.faa ECOR_27_genome.faa ECOR_37_genome.faa REL606.fasta sequence_RS.fasta sequence_oldMG1655.fasta | sed 's/\([ACGT]\)\(>\)/\1\n\2/g' > Ecoli_all_genomes.faa
cat Ecoli_annotated.gtf | grep -v "#" | gtf2bed | bedtools getfasta -s -name -fi Ecoli_all_genomes.faa -bed - > Ecoli_annotated.faa

cat ../Salmonella_list/GCA_000210855.2_ASM21085v2_genomic.fna ../Salmonella_list/Salmonella_14028s.fasta | sed 's/\([ACGT]\)\(>\)/\1\n\2/g' > Salmonella_all_genomes.faa
cat Salmonella_annotated.gtf | grep -v "#" | gtf2bed | bedtools getfasta -s -name -fi Salmonella_all_genomes.faa -bed - > Salmonella_annotated.faa

cat ../Mycobacterium_list/H37Rv.fna > Mycobacterium_all_genomes.faa
cat Mycobacterium_annotated.gtf | grep -v "#" | gtf2bed | bedtools getfasta -s -name -fi Mycobacterium_all_genomes.faa -bed - > Mycobacterium_annotated.faa

#non-red with usearch

for i in Ecoli Salmonella Mycobacterium
do
usearch -sortbylength "$i"_annotated.faa -fastaout "$i"_annotated.sorted.faa -minseqlength 1
usearch -cluster_smallmem "$i"_annotated.sorted.faa -id 0.9 -centroids "$i"_annotated.nr.faa -uc "$i"_clusters.uc
done





#for i in Ecoli Salmonella Mycobacterium
#do
#for j in $(cut -f1 "$i"_annotated.gtf | sort -u)
#do
#awk -F '\t' '($7=="+")' "$i"_annotated.gtf | cut -f5 | sort -u > "$i"_annotated_stops_plus
#awk -F '\t' '($7=="-")' "$i"_annotated.gtf | cut -f4 | sort -u > "$i"_annotated_stops_minus
#done
#done

#for j in $(cat "$i"_annotated_stops_plus); do awk -v var="$j" -F '\t' '($7=="+"&&$5==var)' "$i"_annotated.gtf | sort -nk 4 | head -1; done > "$i"_annotated.nr.gtf
#awk -F '\t' '($7=="-")' "$i"_annotated.gtf | cut -f4 | sort -u > "$i"_annotated_stops_minus
#for j in $(cat "$i"_annotated_stops_minus); do awk -v var="$j" -F '\t' '($7=="-"&&$4==var)' "$i"_annotated.gtf | sort -nrk 5 | head -1; done >> "$i"_annotated.nr.gtf
#done

#cat "$i"_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | seqkit fx2tab | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/^/>/g" | sed "s/\t$//g" | sed "s/\t/\n/g" > "$i"_all_annotated.nr.seqkit.faa
#cat "$i"_annotated.nr.gtf | gtf2bed | bedtools getfasta -s -name -fi "$i".faa -bed - | seqkit rmdup -s - | grep "^>" | tr -d ">" | cut -f1 -d "(" | sed "s/.*/\"&\"/g" | grep -F -f - "$i"_annotated.nr.gtf > "$i"_all_annotated.nr.seqkit.gtf

#Get their sequences and proteins
#cat Ecoli_all_protogenes_final.gtf | gtf2bed | bedtools getfasta -s -name -fi all_K12_variants.faa -bed - > Ecoli_all_protogenes_final.faa
#/stor/work/Ochman/hassan/tools/faTrans -stop Ecoli_all_protogenes_final.faa Ecoli_all_protogenes_final.prot.faa
