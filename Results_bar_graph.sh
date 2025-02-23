#Total novel genes:

cd /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties

cut -f1 -d "(" protogene_exclude.txt | sed "s/.*/\"&\"/g" | grep -v -F -f - Ecoli_queryfile.gtf | egrep -iv "smorf|prod|gms2|balrog" | wc -l
cut -f1 -d "(" protogene_exclude.txt | sed "s/.*/\"&\"/g" | grep -v -F -f - Salmonella_queryfile.gtf | egrep -iv "smorf|prod|gms2|balrog" | wc -l
cut -f1 -d "(" protogene_exclude.txt | sed "s/.*/\"&\"/g" | grep -v -F -f - Mycobacterium_queryfile.gtf | egrep -iv "smorf|prod|gms2|balrog" | wc -l

#Intergenic?
cut -f1 -d "(" protogene_exclude.txt | sed "s/.*/\"&\"/g" | grep -v -F -f - Ecoli_queryfile.gtf | #all genes
egrep -iv "smorf|prod|gms2|balrog" | cut -f2 -d "\"" | sort -u | sed "s/$/(/g" | #all novel genes
grep -F -f - ../Ecoli_CDS_queryfile.mapped.REL606.gtf | #all mappable
bedtools intersect -v -a - -b ../Ecoli_nonredundant_REL606.gtf | wc -l

cut -f1 -d "(" protogene_exclude.txt | sed "s/.*/\"&\"/g" | grep -v -F -f - Salmonella_queryfile.gtf | #all genes
egrep -iv "smorf|prod|gms2|balrog" | cut -f2 -d "\"" | sort -u | sed "s/$/(/g" | #all novel genes
grep -F -f - ../Salmonella.tobemapped.protogenes.gtf | #all mappable
bedtools intersect -v -a - -b ../Salmonella_nonredundant.gtf |
wc -l

cut -f1 -d "(" protogene_exclude.txt | sed "s/.*/\"&\"/g" | grep -v -F -f - Mycobacterium_queryfile.gtf | #all genes
egrep -iv "smorf|prod|gms2|balrog" | #all novel genes
bedtools intersect -v -a - -b ../Mycobacterium_nonredundant_annotated.gtf |
wc -l

#Total annotated genes:

wc -l ../Ecoli_nonredundant_REL606.gtf
wc -l ../Salmonella_nonredundant.gtf
wc -l ../Mycobacterium_nonredundant_annotated.gtf

#Genus-specific novel genes:

egrep -iv "smorf|prod|gms2|balrog" Ecoli_genusspecific_ORFans.final.txt | wc -l
egrep -iv "smorf|prod|gms2|balrog" Salmonella_genusspecific_ORFans.final.txt | wc -l
egrep -iv "smorf|prod|gms2|balrog" Mycobacterium_speciesspecific_ORFans.final.txt | wc -l

#Intergenic?

egrep -iv "smorf|prod|gms2|balrog" Ecoli_genusspecific_ORFans.final.txt | #novel genes
grep -F -f - ../Ecoli_CDS_queryfile.mapped.REL606.gtf | #all mappable
bedtools intersect -v -a - -b ../Ecoli_nonredundant_REL606.gtf | wc -l

egrep -iv "smorf|prod|gms2|balrog" Salmonella_genusspecific_ORFans.final.txt | #novel genes
grep -F -f - ../Salmonella.tobemapped.protogenes.gtf | #all mappable
bedtools intersect -v -a - -b ../Salmonella_nonredundant.gtf | wc -l

egrep -iv "smorf|prod|gms2|balrog" Mycobacterium_speciesspecific_ORFans.final.txt | cut -f1 -d "(" | sed "s/.*/\"&\"/g" |
grep -F -f - Mycobacterium_queryfile.gtf |
bedtools intersect -v -a - -b ../Mycobacterium_nonredundant_annotated.gtf | wc -l

#Genus-specific annotated genes:

egrep -i "smorf|prod|gms2|balrog" Ecoli_genusspecific_ORFans.final.txt | grep -v "NZ_" | wc -l
egrep -i "smorf|prod|gms2|balrog" Salmonella_genusspecific_ORFans.final.txt | wc -l
egrep -i "smorf|prod|gms2|balrog" Mycobacterium_speciesspecific_ORFans.final.txt | wc -l

#Traceable: total, novel
cd /stor/work/Ochman/hassan/protogene_extension/comparative_genomics/flanks
grep -v -w -F -f ../Ecoli_step1_speciespecific_ORFan.txt ../Ecoli_step1_genusspecific_ORFan.txt | cut -f1 -d "(" | sed "s/^/ls Ecoli*/g" | sed "s/$/_compiled.intervalinfo.txt/g" | bash | rev | cut -f2- -d "_" | rev | cut -f2- -d "_" > /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/Ecoli_genusspecific_traceable.txt
egrep -iv "smorf|prod|gms2|balrog" Ecoli_genusspecific_ORFans.final.txt | cut -f1 -d "(" | grep -w -F -f - Ecoli_genusspecific_traceable.txt | wc -l

grep -v -w -F -f ../Salmonella_step1_speciespecific_ORFan.txt ../Salmonella_step1_genusspecific_ORFan.txt | cut -f1 -d "(" | sed "s/^/ls Salmonella*/g" | sed "s/$/_compiled.intervalinfo.txt/g" | bash | rev | cut -f2- -d "_" | rev | cut -f2- -d "_" > /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/Salmonella_genusspecific_traceable.txt
egrep -iv "smorf|prod|gms2|balrog" Salmonella_genusspecific_ORFans.final.txt | cut -f1 -d "(" | grep -w -F -f - Salmonella_genusspecific_traceable.txt | wc -l

ls Myc*_compiled.intervalinfo.txt | rev | cut -f2- -d "_" | rev | cut -f2- -d "_" | sort -u > /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/Mycobacterium_speciesspecific_traceable.txt
egrep -iv "smorf|prod|gms2|balrog" Mycobacterium_speciesspecific_ORFans.final.txt | cut -f1 -d "(" | grep -w -F -f - Mycobacterium_total_traceable.txt | wc -l

cat ../Ecoli_step1_genusspecific_ORFan.txt | cut -f1 -d "(" | sed "s/^/ls Ecoli*/g" | sed "s/$/_compiled.intervalinfo.txt/g" | bash | rev | cut -f2- -d "_" | rev | cut -f2- -d "_" > /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/Ecoli_speciesspecific_traceable.txt


#Species-specific ORFans:
grep -v -F -f ../../comparative_genomics/Ecoli_step1_speciespecific_nonORFan.txt Ecoli_genusspecific_ORFans.final.txt > Ecoli_speciesspecific_ORFans.step2.txt
grep "non-ORFan sp" alignment_info  | cut -f1 | rev | cut -f2- -d "_" | rev | cut -f2- -d "_" | sed "s/$/(/g" | grep -v -F -f - Ecoli_speciesspecific_ORFans.step2.txt > Ecoli_speciesspecific_ORFans.final.txt
cd /stor/work/Ochman/hassan/protogene_extension/comparative_genomics/flanks
cat /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/Ecoli_speciesspecific_ORFans.final.txt | egrep -iv "prod|smorf|balrog|gms2" | cut -f1 -d "(" | sed "s/^/ls Ecoli_/g" | sed "s/$/_compiled.intervalinfo.txt/g" | bash

grep -v -F -f ../../comparative_genomics/Salmonella_step1_speciespecific_nonORFan.txt Salmonella_genusspecific_ORFans.final.txt > Salmonella_speciesspecific_ORFans.step2.txt
grep "non-ORFan sp" alignment_info | cut -f1 | rev | cut -f2- -d "_" | rev | cut -f2- -d "_" | sed "s/$/(/g" | grep -v -F -f - Salmonella_speciesspecific_ORFans.step2.txt > Salmonella_speciesspecific_ORFans.final.txt
cd /stor/work/Ochman/hassan/protogene_extension/comparative_genomics/flanks
cat /stor/work/Ochman/hassan/protogene_extension/expression_location_properties/sequence_properties/Salmonella_speciesspecific_ORFans.final.txt | egrep -iv "prod|smorf|balrog|gms2" | cut -f1 -d "(" | sed "s/^/ls Salmonella_/g" | sed "s/$/_compiled.intervalinfo.txt/g" | bash


