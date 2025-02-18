mkdir /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties
cd /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties

cp /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/*query*faa .
cp /stor/work/Ochman/hassan/protogene_extension/comparative_genomics/*query*gtf .

#Get control ORFs for comparison
cp /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/sequence_RS.fasta Ecoli.focalstrain.faa 
cp /stor/work/Ochman/hassan/protogene_extension/Salmonella_list/GCA_000210855.2_ASM21085v2_genomic.fna Salmonella.focalstrain.faa
cp /stor/work/Ochman/hassan/protogene_extension/Mycobacterium_list/H37Rv.fna Mybacterium.focalstrain.faa


for i in Ecoli Salmonella Mycobacterium
do
getorf -sequence "$i".focalstrain.faa -outseq "$i".getorf.bacterial -table 1 -minsize 30 -find 3
seqkit fx2tab "$i".getorf.bacterial | sed "s/\t$//g" | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/g" | grep "^>" | grep -v "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$2,"\t",$4+3,"\t.\t+\t0\ttranscript_id \"",$1,"\";gene_id \"",$1,"\";"}' | sed "s/_/\t/1" | cut -f1,3- > "$i".getorf.bacterial.gtf
seqkit fx2tab "$i".getorf.bacterial | sed "s/\t$//g" | grep -P -v "\tCTG" | sed "s/^/>/g" | sed "s/\t/\n/g" | grep "^>" | grep "REVERSE" | sed "s/\[//g" | sed "s/\]//g" | tr -d ">" | awk '{OFS=""}{print $1"\t.\tCDS\t",$4-3,"\t",$2,"\t.\t-\t0\ttranscript_id \"",$1,"\";gene_id \"",$1,"\";"}' | sed "s/_/\t/1" | cut -f1,3- >> "$i".getorf.bacterial.gtf
done

sed -i "s/^NC/NC_999999.9/g" Mycobacterium.getorf.bacterial.gtf
sed -i "s/^NC/NC_000913.3/g" Ecoli.getorf.bacterial.gtf

#Annotated gtfs
for i in Ecoli Salmonella Mycobacterium
do
cut -f1 "$i".getorf.bacterial.gtf | sort -u | grep -F -f - /stor/work/Ochman/hassan/protogene_extension/Ecoli_list/"$i"_annotated.gtf | awk '
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
}' > "$i".annotated.gtf
done

for i in Ecoli Samonella Mycobacterium
do
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap_build -D . -d "$i" "$i".focalstrain.faa
done

#Ecoli:
#The ones already extracted from focal strain doesn't need to be mapped:
egrep -iv "prodigal|smorf|gms2|balrog" Ecoli_queryfile.gtf | grep -v "^NZ_" | egrep "^NC_000913.3" > Ecoli.mapped.protogenes.gtf
#For others:
egrep -iv "prodigal|smorf|gms2|balrog" Ecoli_queryfile.gtf | grep -v "^NZ_" | egrep -v "^NC_000913.3" | cut -f2 -d "\"" | sed "s/$/(/g" | grep --no-group-separator -A1 -F -f - Ecoli_CDS_queryfile.faa > Ecoli.tobemapped.protogenes.faa
/stor/work/Ochman/hassan/tools/gmap-2021-05-27/bin/gmap -D . -d Ecoli -f 2 --gff3-fasta-annotation=1 Ecoli.tobemapped.protogenes.faa > Ecoli.tobemapped.protogenes.gff3
awk -F '\t' '($3=="mRNA")' Ecoli.tobemapped.protogenes.gff3 | grep "mrna1" | cut -f1 -d ";" | awk -F '\t' '{OFS=FS}{print $1,$2,"CDS",$4,$5,$6,$7,$8,"transcript_id \""$9"\";gene_id \""$9"\";"}' | sed "s/.mrna1//g" | sed "s/ID=//g" > Ecoli.tobemapped.protogenes.gtf

#Get controls:
cat Ecoli.annotated.gtf Ecoli.mapped.protogenes.gtf Ecoli.tobemapped.protogenes.gtf | bedtools intersect -s -v -a Ecoli.getorf.bacterial.gtf -b - | gtf2bed | bedtools getfasta -s -name -fi Ecoli.focalstrain.faa -bed - | sed "s/>/>control_/g" > Ecoli.control.CDS.faa
egrep -iv "prodigal|smorf|gms2|balrog" Salmonella_queryfile.gtf | cat Salmonella.annotated.gtf - | bedtools intersect -s -v -a Salmonella.getorf.bacterial.gtf -b - | gtf2bed | bedtools getfasta -s -name -fi Salmonella.focalstrain.faa -bed - | sed "s/>/>control_/g" > Salmonella.control.CDS.faa
egrep -iv "prodigal|smorf|gms2|balrog" Mybacterium_queryfile.gtf | cat Mycobacterium.annotated.gtf - | bedtools intersect -s -v -a Mycobacterium.getorf.bacterial.gtf -b - | gtf2bed | bedtools getfasta -s -name -fi Mycobacterium.focalstrain.faa -bed - | sed "s/>/>control_/g" > Mycobacterium.control.CDS.faa

for i in Ecoli Salmonella Mycobacterium; do /stor/work/Ochman/hassan/tools/faTrans -stop "$i".control.CDS.faa "$i".control.prot.faa; done

for i in Ecoli Salmonella Mycobacterium
do
cat "$i"_CDS_queryfile.faa "$i".control.CDS.faa > "$i"_CDS.final.faa
cat "$i"_protein_queryfile.faa "$i".control.prot.faa > "$i"_prot.final.faa
done

for i in Ecoli Salmonella Mycobacterium
do
#length
cat "$i"_CDS.final.faa | seqkit fx2tab | awk -F '\t' '{OFS=","}{print $1,length($2)}' | tr -d ">" | sed "s/$/,length/g" > "$i"_sequenceproperties.csv
#GC
cat "$i"_CDS.final.faa | seqkit fx2tab | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" | awk '$0 ~ ">" { if (NR > 1) { print name, gc_count / total_length * 100 "%"; } name = substr($0, 2); total_length = 0; gc_count = 0; next; } { total_length += length($0); gc_count += gsub(/[GCgc]/, "", $0); } END { print name, gc_count / total_length * 100 "%"; }' | sed "s/%//g" | sed "s/ /,/g" | sed "s/,,/,/g" | sed "s/$/,GC/g" >> "$i"_sequenceproperties.csv
#GC_3rd
cat "$i"_CDS.final.faa | seqkit fx2tab | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" | awk '$0 ~ ">" { if (NR > 1) { print name, third_gc_count / (total_length / 3) * 100 "%"; } name = substr($0, 2); total_length = 0; third_gc_count = 0; next; } { total_length += length($0); for (i = 3; i <= length($0); i += 3) { if (substr($0, i, 1) ~ /[GCgc]/) { third_gc_count++; } } } END { print name, third_gc_count / (total_length / 3) * 100 "%"; }' | sed "s/%//g" | sed "s/ /,/g" | sed "s/,,/,/g" | sed "s/$/,GC_3rd/g" >> "$i"_sequenceproperties.csv
#Polar AA content
cat "$i"_CDS.final.faa | seqkit fx2tab | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" | awk '$0 ~ ">" { if (NR > 1) { print name, polar_count / total_length * 100 "%"; } name = substr($0, 2); total_length = 0; polar_count = 0; next; } { total_length += length($0); polar_count += gsub(/[STYECQNHKR]/, "", $0); } END { print name, polar_count / total_length * 100 "%"; }' | sed "s/%//g" | sed "s/ /,/g" | sed "s/,,/,/g" | sed "s/$/,polarAA/g" >> "$i"_sequenceproperties.csv
#Hydrophobic AA content
cat "$i"_CDS.final.faa | seqkit fx2tab | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" | awk '$0 ~ ">" { if (NR > 1) { print name, polar_count / total_length * 100 "%"; } name = substr($0, 2); total_length = 0; polar_count = 0; next; } { total_length += length($0); polar_count += gsub(/[ACFGILMPVWY]/, "", $0); } END { print name, polar_count / total_length * 100 "%"; }' | sed "s/%//g" | sed "s/ /,/g" | sed "s/,,/,/g" | sed "s/$/,hydrophobicAA/g" >> "$i"_sequenceproperties.csv
#Biosynthetic cost - Akashi and Gojobari lists
cat "$i"_CDS.final.faa | seqkit fx2tab | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" | awk 'NR==FNR{costs[$1]=$2; next} $0 ~ ">" { if (NR > 1 && total_length > 0) { print name, total_cost / total_length; } name = substr($0, 2); total_length = 0; total_cost = 0; next; } { for (i = 1; i <= length($0); i++) { if (substr($0, i, 1) in costs) { total_cost += costs[substr($0, i, 1)]; } } total_length += length($0); } END { if (total_length > 0) { print name, total_cost / total_length; } }' /stor/work/Ochman/hassan/tools/AA_metabolic_costs.txt - | sed "s/%//g" | sed "s/ /,/g" | sed "s/,,/,/g" | sed "s/$/,metabol/g" >> "$i"_sequenceproperties.csv
done
#cai
/stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/cai -seqall Ecoli_CDS.final.faa -cfile /stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/data/CODONS/Eecoli.cut -outfile Ecoli.cai
/stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/cai -seqall Salmonella_CDS.final.faa -cfile /stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/data/CODONS/Esalty.cut -outfile Salmonella.cai
/stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/cai -seqall Mycobacterium_CDS.final.faa -cfile /stor/work/Ochman/hassan/tools/EMBOSS-6.6.0/emboss/data/CODONS/Emyctu.cut -outfile Mycobacterium.cai
for i in Ecoli Salmonella Mycobacterium
do
cut -f2,4 -d " " "$i".cai | sed "s/ /,/g" | sed "s/$/,cai/g" >> "$i"_sequenceproperties.csv
done

#Assign datasets

egrep -i "balrog|smorf|prodigal|gms2" Ecoli_sequenceproperties.csv | grep -v "^NZ_" | sed "s/$/,annotated,annotated/g" > Ecoli_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Ecoli_sequenceproperties.csv | grep -v "^NZ_" | grep "Stringer" | sed "s/$/,novel,Stringer/g" >> Ecoli_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Ecoli_sequenceproperties.csv | grep -v "^NZ_" | grep "Nakahigashi" | sed "s/$/,novel,Nakahigashi/g" >> Ecoli_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Ecoli_sequenceproperties.csv | grep -v "^NZ_" | grep "VanOrsdel" | sed "s/$/,novel,VanOrsdel/g" >> Ecoli_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Ecoli_sequenceproperties.csv | grep -v "^NZ_" | grep "Weaver" | sed "s/$/,novel,Weaver/g" >> Ecoli_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Ecoli_sequenceproperties.csv | grep -v "^NZ_" | grep "^EC" | sed "s/$/,novel,Ndah/g" >> Ecoli_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Ecoli_sequenceproperties.csv | grep -v "^NZ_" | grep "^NC" | sed "s/$/,novel,MS/g" >> Ecoli_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Ecoli_sequenceproperties.csv | grep -v "^NZ_" | grep "control" | sed "s/$/,control,control/g" >> Ecoli_sequenceproperties.marked.csv

egrep -i "balrog|smorf|prodigal|gms2" Salmonella_sequenceproperties.csv | sed "s/$/,annotated,annotated/g" > Salmonella_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Salmonella_sequenceproperties.csv | grep "Fijalkowski" | sed "s/$/,novel,Fijalkowski/g" >> Salmonella_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Salmonella_sequenceproperties.csv | grep "Giess" | sed "s/$/,novel,Giess/g" >> Salmonella_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Salmonella_sequenceproperties.csv | grep "Ndah" | sed "s/$/,novel,Ndah/g" >> Salmonella_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Salmonella_sequenceproperties.csv | grep "Venturini" | sed "s/$/,novel,Venturini/g" >> Salmonella_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Salmonella_sequenceproperties.csv | grep "control" | sed "s/$/,control,control/g" >> Salmonella_sequenceproperties.marked.csv

egrep -i "balrog|smorf|prodigal|gms2" Mycobacterium_sequenceproperties.csv | sed "s/$/,annotated,annotated/g" > Mycobacterium_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Mycobacterium_sequenceproperties.csv | grep -i "Riboret" | sed "s/$/,novel,RiboRET/g" >> Mycobacterium_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Mycobacterium_sequenceproperties.csv | grep -i "leaderless" | sed "s/$/,novel,Leaderless/g" >> Mycobacterium_sequenceproperties.marked.csv
egrep -iv "balrog|smorf|prodigal|gms2" Mycobacterium_sequenceproperties.csv | grep "control" | sed "s/$/,control,control/g" >> Mycobacterium_sequenceproperties.marked.csv


#Now categorize
#Get the ORFans first
sed "s/\t m/\tnon-ORFan/g" alignment_info | sed "s/\tm/\tnon-ORFan/g" | sed "s/_mafft.aln//g" | awk -F '\t' '($2=="non-ORFan")' | grep "^Ecoli" | cut -f1 | sort -u | sed "s/Ecoli_//g" | sed "s/$/(/g" | grep -F -vf - ../../comparative_genomics/Ecoli_step1_genusspecific_ORFan.txt | grep -v "^NZ_" > Ecoli_genusspecific_ORFans.txt
cat ../../comparative_genomics/*_vs_pangenome_annotated.tsv | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1 | egrep -iv "NZ_|prodigal|gms2|balrog|smorf" | sort -u > protogene_exclude.txt
grep -v -F -f protogene_exclude.txt Ecoli_genusspecific_ORFans.txt > Ecoli_genusspecific_ORFans.final.txt

sed "s/\t m/\tnon-ORFan/g" alignment_info | sed "s/\tm/\tnon-ORFan/g" | sed "s/_mafft.aln//g" | awk -F '\t' '($2=="non-ORFan")' | grep "^Salmonella" | cut -f1 | sort -u | sed "s/Salmonella_//g" | sed "s/$/(/g" | grep -F -vf - ../../comparative_genomics/Salmonella_step1_genusspecific_ORFan.txt > Salmonella_genusspecific_ORFans.txt
grep -v -F -f protogene_exclude.txt Salmonella_genusspecific_ORFans.txt > Salmonella_genusspecific_ORFans.final.txt

sed "s/\t m/\tnon-ORFan/g" alignment_info | sed "s/\tm/\tnon-ORFan/g" | sed "s/_mafft.aln//g" | awk -F '\t' '($2=="non-ORFan")' | grep "^Mycobacterium" | cut -f1 | sort -u | sed "s/Mycobacterium_//g" | sed "s/$/(/g" | grep -F -vf - ../../comparative_genomics/Mycobacterium_step1_speciespecific_ORFan.txt > Mycobcaterium_speciesspecific_ORFans.txt
grep -v -F -f protogene_exclude.txt Mycobcaterium_speciesspecific_ORFans.txt > Mycobcaterium_speciesspecific_ORFans.final.txt

cat genomes/accession_genomeID_taxonomy.tsv proteins/accession_proteinID_taxonomy.tsv | grep -v -F -f odd_taxonomic_names.txt - | sed "s/\[//g" | sed "s/\t\'/\t/g" | sed "s/\tuncultured_/\t/g" | sed "s/\tCandidatus_/\t/g" | grep -P -v "\tbacterium" | grep -P -v "\tcandidate_" | cut -f-2 -d "_" | cut -f2- > GBRS_ids_taxa.tsv
sort -k1 /stor/scratch/Ochman/hassan/100724_Complete_Genomes/GBRS_ids_taxa.tsv -o /stor/scratch/Ochman/hassan/100724_Complete_Genomes/GBRS_ids_taxa.tsv

#Conserved, non-ORFans

cat ../../comparative_genomics/Ecoli_vs_GBRS_annotated.tsv ../../comparative_genomics/Ecoli_vs_GBRS_ORFs.tsv | awk -F '\t' '($5>60&&$16<0.001)' | cut -f1,2 | sort -u > Ecoli_extragenus_hits.tsv

cat ../../comparative_genomics/Salmonella_vs_GBRS_annotated.tsv | awk -F '\t' '($16<0.001&&$5>60)' | cut -f1,2 > Sal_test1
cat ../../comparative_genomics/Salmonella_vs_GBRS_ORFs.tsv | awk -F '\t' '($16<0.001&&$5>60)' | cut -f1,2 > Sal_test2
rev Sal_test2 | cut -f2- -d "_" | rev | cat Sal_test1 - > Salmonella_extragenus_hits.tsv
sort -u Salmonella_extragenus_hits.tsv -o Salmonella_extragenus_hits.tsv

cat ../../comparative_genomics/Mycobacterium_vs_GBRS_annotated.tsv | awk -F '\t' '($16<0.001&&$5>60)' | cut -f1,2 > Myc_test1
cat ../../comparative_genomics/Mycobacterium_vs_GBRS_ORFs.tsv | awk -F '\t' '($16<0.001&&$5>60)' | cut -f1,2 > Myc_test2
rev Myc_test2 | cut -f2- -d "_" | rev | cat Myc_test1 - > Mycobacterium_extragenus_hits.tsv
sort -u Mycobacterium_extragenus_hits.tsv -o Mycobacterium_extragenus_hits.tsv

sed -i "s/^/\t/g" odd_taxonomic_names.txt

awk '{sub(/_.*/, "", $2)}1' Ecoli_extragenus_hits.tsv | sort -k2 | join -1 2 -2 1 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/GBRS_ids_taxa.tsv > Ecoli_extragenus_hits.taxa.tsv
###
awk '{sub(/_.*/, "", $2)}1' Ecoli_extragenus_hits.tsv > Ecoli_interim
###
sort -k2 Salmonella_extragenus_hits.tsv | join -1 2 -2 1 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/GBRS_ids_taxa.tsv > Salmonella_extragenus_hits.taxa.tsv
sort -k2 Mycobacterium_extragenus_hits.tsv | join -1 2 -2 1 - /stor/scratch/Ochman/hassan/100724_Complete_Genomes/GBRS_ids_taxa.tsv > Mycobacterium_extragenus_hits.taxa.tsv

cut -f2- -d " " Ecoli_extragenus_hits.taxa.tsv | sort -u | cut -f1 -d " " | sort | uniq -c > Ecoli_extragenus_hits.genusnumber
cut -f2- -d " " Salmonella_extragenus_hits.taxa.tsv | sort -u | cut -f1 -d " " | sort | uniq -c > Salmonella_extragenus_hits.genusnumber
cut -f2- -d " " Mycobacterium_extragenus_hits.taxa.tsv | sort -u | cut -f1 -d " " | sort | uniq -c > Mycobacterium_extragenus_hits.genusnumber

for i in Ecoli Salmonella Mycobacterium
do
annot_median=$(egrep "balrog|prodigal|smorfer|gms2" "$i"_extragenus_hits.genusnumber | awk '{print $1}' | sort -n | awk ' { a[i++] = $1; } END { print (i % 2 == 1) ? a[int(i/2)] : (a[int(i/2)-1] + a[int(i/2)]) / 2; }')
novel_median=$(egrep -iv "balrog|prodigal|smorfer|gms2" "$i"_extragenus_hits.genusnumber | awk '{print $1}' | sort -n | awk ' { a[i++] = $1; } END { print (i % 2 == 1) ? a[int(i/2)] : (a[int(i/2)-1] + a[int(i/2)]) / 2; }')
egrep "balrog|prodigal|smorfer|gms2" "$i"_extragenus_hits.genusnumber | awk -v var=$annot_median '($1>var)' | rev | cut -f1 -d " " | rev | sort -u > "$i"_annot_conserved
egrep -iv "balrog|prodigal|smorfer|gms2" "$i"_extragenus_hits.genusnumber | awk -v var=$novel_median '($1>var)' | rev | cut -f1 -d " " | rev | sort -u > "$i"_novel_conserved
done



