cd /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/supp_allprotein_categories

#Files:

Mycobacterium_602_protogenes.prot.sorted.faa  Mycobacterium_602_protogenes.CDS.faa  Salmonella_135_protogenes.prot.sorted.faa  Salmonella_135_protogenes.CDS.faa
Mycobacterium_602_protogenes.prot.faa         Ecoli_559_protogenes.prot.sorted.faa  Salmonella_135_protogenes.prot.faa         Ecoli_559_protogenes.CDS.faa

for i in Ecoli Salmonella Mycobacterium
do
#usearch -cluster_smallmem "$i"_*_protogenes.prot.sorted.faa -id 0.9 -centroids "$i"_protogenes.prot.sorted.nr.faa -uc "$i"_protogenes.prot.clusters.uc
seqkit fx2tab "$i"_protogenes.prot.sorted.nr.faa | sed "s/\t$//g" | sed "s/^/>/g" | sed "s/\t/\n/g" > temp && mv temp "$i"_protogenes.prot.sorted.nr.faa
done

#Load up info on dataset, ID, species, strain, method
grep "^>" Ecoli_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tE. coli/g" | egrep -i "Stringer|Weaver|Nakahigashi|EC_Ndah" | awk -F "_" '{print $1"\t"$0}' | sed "s/^EC/EC_Ndah/g" | sed "s/$/\tK-12 MG1655/g" | sed "s/$/\tRiboSeq/g" > master_table_interim.tsv
grep "^>" Ecoli_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tE. coli/g" | egrep -i "VanOrsdel" | awk -F "_" '{print $1"\t"$0}' | sed "s/$/\tK-12 MG1655/g" | sed "s/$/\tWestern blotting/g" >> master_table_interim.tsv
grep "^>" Ecoli_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tE. coli/g" | egrep -i "NC_012967.1" | sed "s/^/Caglar (re-analysis)\t/g" | sed "s/$/\tREL606/g" | sed "s/$/\tMass spectrometry/g" >> master_table_interim.tsv
grep "^>" Ecoli_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tE. coli/g" | egrep -i "NC_000913.3" | sed "s/^/Mori (re-analysis)\t/g" | sed "s/$/\tK-12 MG1655/g" | sed "s/$/\tMass spectrometry/g" >> master_table_interim.tsv
grep "^>" Ecoli_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tE. coli/g" | egrep -i "NZ_QOWW01000007.1" | sed "s/$/\tECOR_11/g" | sed "s/^/This study\t/g" | sed "s/$/\tMass spectrometry/g" >> master_table_interim.tsv
grep "^>" Ecoli_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tE. coli/g" | egrep -i "NZ_QOXM01000014.1|NZ_QOXM01000017.1|NZ_QOXM01000009.1" | sed "s/$/\tECOR_27/g" | sed "s/^/This study\t/g" | sed "s/$/\tMass spectrometry/g" >> master_table_interim.tsv
grep "^>" Ecoli_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tE. coli/g" | egrep -i "NZ_QOXW01000021.1|NZ_QOXW01000045.1|NZ_QOXW01000027.1|NZ_QOXW01000022.1|NZ_QOXW01000003.1" | sed "s/$/\tECOR_37/g" | sed "s/^/This study\t/g" | sed "s/$/\tMass spectrometry/g" >> master_table_interim.tsv
#OK if this wasn't bad enough, I found out I didn't include Baek for Salmonella
grep "^>" Salmonella_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tS. enterica/g" | awk -F "_" '{print $1"\t"$0}' | sed "s/^ST/ST_Ndah/g" | sed "s/$/\tSL1344/g" | sed "s/$/\tRiboSeq/g" >> master_table_interim.tsv
grep "^>" Mycobacterium_protogenes.prot.sorted.nr.faa | tr -d ">" | cut -f1 -d "(" | sed "s/$/\tM. tuberculosis/g" | awk -F "_" '{print $1"\t"$0}' | sed "s/$/\tH37Rv/g" | sed "s/$/\tRiboSeq/g" >> master_table_interim.tsv

#Watch the tab
sort -t '  ' -k2 master_table_interim.tsv -o master_table_interim.tsv
#Watch the tab
cat master_table_interim.tsv | cut -f2 | sed "s/$/(/g" | grep -A1 --no-group-separator -F -f - *prot.sorted.nr.faa | sed '/^>/! s/-/:/' | cut -f2- -d ":" | seqkit fx2tab | sed "s/(+)//g" | sed "s/(:)//g" | awk -F '\t' '{OFS=FS}{print $1,$2,length($2)}' | sort -k1 | join -t '  ' -1 1 -2 2 - master_table_interim.tsv > temp && mv temp master_table_interim.tsv

#Chromosomes and coordinates:
sort -k1 master_table_interim.tsv -o master_table_interim.tsv
cut -f1 master_table_interim.tsv | sed "s/.*/\"&\"/g" | grep -F -f - ../../../Ecoli_list/all_protogenes.gtf | cut -f1,4,5,7,9 | cut -f-2 -d "\"" | sed "s/transcript_id \"//g" | sed "s/ /\t/g" | sort -k5 | join -t ' ' -1 5 -2 1 - master_table_interim.tsv > temp && sort -k1 temp -o master_table_interim.tsv

#Other datasets:
cut -f1 master_table_interim.tsv | sort -u | grep -F -f - *uc | grep ":H" | rev | cut -f1,2 | rev | egrep -iv "Riboret|leaderless" > multiple_datasets.interim
#Manually remove same dataset hits
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - multiple_datasets.interim | awk '{print $2,$1}' | rev | cut -f2- -d "_" | rev | sed "s/ /\t/g" > multiple_datasets.2.interim
awk -F '\t' '{print $2"\t"$1}' multiple_datasets.2.interim | awk -F'\t' '{ values[$2] = (values[$2] == "" ? $1 : values[$2] ", " $1) } END { for (value in values) { print value "\t" values[value] } }' | sed "s/ //g" | sed "s/(+)//g" | sed "s/(-)//g" | sort -k1 | join -t '  ' -1 1 -2 1 - master_table_interim.tsv > temp
cut -f1 temp | grep -v -w -F -f - master_table_interim.tsv | awk -F '\t' '{OFS=FS}{print $1,"N/A",substr($0, index($0, $2))}' >> temp
mv temp master_table_interim.tsv

cut -f1 master_table_interim.tsv | sed "s/.*/\"&\"/g" | grep -F -f - ../../../Ecoli_list/all_protogenes.gtf | bedtools intersect -v -a - -b ../../../Ecoli_list/annotated.gtf | cut -f2 -d "\"" | sort -u | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tintergenic/g" > temp
cut -f1 master_table_interim.tsv | sed "s/.*/\"&\"/g" | grep -F -f - ../../../Ecoli_list/all_protogenes.gtf | bedtools intersect -v -a - -b ../../../Ecoli_list/annotated.gtf | cut -f2 -d "\"" | sort -u | grep -v -w -F -f - master_table_interim.tsv | sed "s/$/\tnot_intergenic/g" >> temp
mv temp master_table_interim.tsv

#genus specific ORFans:
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Ecoli_genusspecific_ORFans.txt | cut -f1 -d "(" | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tgenus_specific/g" > temp
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Ecoli_genusspecific_ORFans.txt | cut -f1 -d "(" | grep -v -w -F -f - master_table_interim.tsv | grep "E. coli" | sed "s/$/\tgenus_nonORFan/g" >> temp
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Salmonella_genusspecific_ORFans.txt | cut -f1 -d "(" | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tgenus_specific/g" >> temp
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Salmonella_genusspecific_ORFans.txt | cut -f1 -d "(" | grep -v -w -F -f - master_table_interim.tsv | grep "S. enterica" | sed "s/$/\tgenus_nonORFan/g" >> temp
egrep -i "RiboRET|Leader" master_table_interim.tsv | sed "s/$/\tN\/A/g" >> temp
mv temp master_table_interim.tsv

#species specific ORFans:
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Ecoli_speciesspecific_ORFans.final.txt | cut -f1 -d "(" | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tspecies_specific/g" > temp
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Ecoli_speciesspecific_ORFans.final.txt | cut -f1 -d "(" | grep -v -w -F -f - master_table_interim.tsv | grep "E. coli" | sed "s/$/\tspecies_nonORFan/g" >> temp
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Salmonella_speciesspecific_ORFans.final.txt | cut -f1 -d "(" | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tspecies_specific/g" >> temp
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Salmonella_speciesspecific_ORFans.final.txt | cut -f1 -d "(" | grep -v -w -F -f - master_table_interim.tsv | grep "S. enterica" | sed "s/$/\tspecies_nonORFan/g" >> temp
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Mycobacterium_speciesspecific_ORFans.final.txt | cut -f1 -d "(" | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tspecies_specific/g" >> temp
cut -f1 master_table_interim.tsv | sed "s/$/(/g" | grep -F -f - ../Mycobacterium_speciesspecific_ORFans.final.txt | cut -f1 -d "(" | grep -v -w -F -f - master_table_interim.tsv | grep "M. tuberculosis" | sed "s/$/\tspecies_nonORFan/g" >> temp
mv temp master_table_interim.tsv

#Genus specific traceable:
grep "E. coli" master_table_interim.tsv | grep "genus_specific" | cut -f1 | grep -w -F -f - ../Ecoli_genusspecific_traceable.txt | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tgenus_traceable/g" > temp
grep "E. coli" master_table_interim.tsv | grep "genus_specific" | cut -f1 | grep -w -F -f - ../Ecoli_genusspecific_traceable.txt | grep -v -w -F -f - master_table_interim.tsv | grep "E. coli" | grep "genus_specific" | sed "s/$/\tgenus_nontraceable/g" >> temp
grep "S. enterica" master_table_interim.tsv | grep "genus_specific" | cut -f1 | grep -w -F -f - ../Salmonella_genusspecific_traceable.txt | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tgenus_traceable/g" >> temp
grep "S. enterica" master_table_interim.tsv | grep "genus_specific" | cut -f1 | grep -w -F -f - ../Salmonella_genusspecific_traceable.txt | grep -v -w -F -f - master_table_interim.tsv | grep "S. enterica" | grep "genus_specific" | sed "s/$/\tgenus_nontraceable/g" >> temp
cut -f1 temp | grep -v -w -F -f - master_table_interim.tsv |  sed "s/$/\tN\/A/g" >> temp
mv temp master_table_interim.tsv

#Species specific traceable:
grep "E. coli" master_table_interim.tsv | grep "species_specific" | cut -f1 | grep -w -F -f - ../Ecoli_speciesspecific_traceable.txt | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tspecies_traceable/g" > temp
grep "E. coli" master_table_interim.tsv | grep "species_specific" | cut -f1 | grep -w -F -f - ../Ecoli_speciesspecific_traceable.txt | grep -v -w -F -f - master_table_interim.tsv | grep "E. coli" | grep "species_specific" | sed "s/$/\tspecies_nontraceable/g" >> temp
grep "M. tuberculosis" master_table_interim.tsv | grep "species_specific" | cut -f1 | grep -w -F -f - ../Mycobacterium_speciesspecific_traceable.txt | grep -w -F -f - master_table_interim.tsv | sed "s/$/\tspecies_traceable/g" >> temp
grep "M. tuberculosis" master_table_interim.tsv | grep "species_specific" | cut -f1 | grep -w -F -f - ../Mycobacterium_speciesspecific_traceable.txt | grep -v -w -F -f - master_table_interim.tsv | grep "tuberculosis" | grep "species_specific" | sed "s/$/\tspecies_nontraceable/g" >> temp
cut -f1 temp | grep -v -w -F -f - master_table_interim.tsv |  sed "s/$/\tN\/A/g" >> temp
mv temp master_table_interim.tsv

#Extragenus hit numbers:
grep "E. coli" master_table_interim.tsv | cut -f1 | sed "s/$/(/g" | grep -F -f - ../Ecoli_extragenus_hits.genusnumber | awk '{print $2,$1}' | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/ /\t/g" > conservation_interim
grep "S. enterica" master_table_interim.tsv | cut -f1 | sed "s/$/(/g" | grep -F -f - ../Salmonella_extragenus_hits.genusnumber | awk '{print $2,$1}' | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/ /\t/g" >> conservation_interim
grep "M. tuberculosis" master_table_interim.tsv | cut -f1 | sed "s/$/(/g" | grep -F -f - ../Mycobacterium_extragenus_hits.genusnumber | awk '{print $2,$1}' | sed "s/(+)//g" | sed "s/(-)//g" | sed "s/ /\t/g" >> conservation_interim
sort -k1 conservation_interim -o conservation_interim
sort -k1 master_table_interim.tsv | join -t '  ' -1 1 -2 1 - conservation_interim
cut -f1 conservation_interim | grep -v -w -F -f - master_table_interim.tsv | sed "s/$/\tN\/A/g" >> temp
mv temp master_table_interim.tsv




