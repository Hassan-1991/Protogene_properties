cd /stor/scratch/Ochman/hassan/112724_protogene_extension

awk '{OFS=""}{print "time java -Xmx3500M -jar \/stor\/work\/Ochman\/hassan\/Fall_2022\/massspec_database\/MSGFplus\/MSGFPlus.jar -s \/stor\/work\/Ochman\/hassan\/mass_spec\/raw_data\/",$1," -d \/stor\/work\/Ochman\/hassan\/proteomics_denovo\/1031_RNAseq_databasemaking\/",$1,"_htseq.tsv.database.tpm1.proteins.fa -inst 1 -t 10ppm -ti 0,1 -mod \/stor\/work\/Ochman\/hassan\/Fall_2022\/08112023_MSGFpercolator\/mods.txt -ntt 2 -tda 1 -maxMissedCleavages 2 -addFeatures 1 -thread 72"}' MURIfiles_with_RNAprotein_data.txt | sed "s/\.mgf//g" > MSGF.sh

#Convert to tsv
cd /stor/work/Ochman/hassan/Fall_2022/massspec_database/MSGFplus/ #For some reason this code only runs in the directory of the jar file
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI_*/*mzid | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i "$1" -o "$1".tsv -showQValue 1 -showDecoy 1 -unroll 1"}' | sed "s/mzid.tsv/tsv/g" > mzid2tsv.sh


ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI_*/*mzid | rev | cut -f 1 -d '/' | cut -f 2- -d '.' | rev | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i /stor/work/Ochman/hassan/proteomics_denovo/htseq/renew/",$1," -o /stor/work/Ochman/hassan/proteomics_denovo/htseq/renew/",$1,".tsv -showQValue 1 -showDecoy 1 -unroll 1"}' > mzidtotsv.sh

ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI_*/*mzid | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i /stor/work/Ochman/hassan/proteomics_denovo/htseq/renew/",$1," -o /stor/work/Ochman/hassan/proteomics_denovo/htseq/renew/",$1,".tsv -showQValue 1 -showDecoy 1 -unroll 1"}'
#ms2rescore

#Mori2021:
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/*mgf | awk '{print "time java -Xmx3500M -jar /stor/work/Ochman/hassan/Fall_2022/massspec_database/MSGFplus/MSGFPlus.jar -s "$0" -d /stor/scratch/Ochman/hassan/112724_protogene_extension/K12MG1655.final.prot.faa -inst 1 -t 10ppm -ti 0,1 -mod /stor/work/Ochman/hassan/Fall_2022/08112023_MSGFpercolator/mods.txt -ntt 2 -tda 1 -maxMissedCleavages 2 -addFeatures 1 -thread 72"}' > running.sh
cd /stor/work/Ochman/hassan/Fall_2022/massspec_database/MSGFplus/ #For some reason this code only runs in the directory of the jar file
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/*mzid | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i "$1" -o "$1".tsv -showQValue 1 -showDecoy 1 -unroll 1"}' | sed "s/mzid.tsv/tsv/g" > running.sh

for i in $(ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/*mzid | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
do
echo "python3 /stor/work/Ochman/hassan/mass_spec/test_files/mgf_search_result_annotator_test1.py --format MSGF_ident --input /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/${i}.mgf --search /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/${i}.mzid --output /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/${i}_annotated.mgf"
done

#Overall
cat *tsv | egrep -i "gms|balrog|smorf|prodigal" | grep -v "XXX" > all_canonical_PSMs.tsv
cat *tsv | egrep "XXX" | cut -f10 | cut -f2- -d '.' | rev | cut -f2- -d '.' | rev | sort -u | grep -o -F -f - all_canonical_PSMs.tsv | sort -u | grep -F -f - *tsv | egrep "XXX" | cut -f11 | sort -u > exclude
cat *tsv | egrep -iv "XXX|gms|balrog|smorf|prodigal" | cut -f10 | cut -f2- -d '.' | rev | cut -f2- -d '.' | rev | sort -u | grep -o -F -f - all_canonical_PSMs.tsv | sort -u | grep -F -f - *tsv | egrep -iv "XXX|gms|balrog|smorf|prodigal" | cut -f11 | sort -u >> exclude

#number of PSMs
wc -l chludwig_Y150505_018_IDA_Lib18-Zhonnge_Ecoli_Lib_18.tsv | cut -f1 -d " "
#number of PSMs passing different cutoffs
for i in 0.01 0.001 0.0001 0.00001
awk -v var="$i" -F '\t' '($16<var)' chludwig_Y150505_018_IDA_Lib18-Zhonnge_Ecoli_Lib_18.tsv | cut -f11 | sort -u | wc -l
done
#Annotated
for i in 0.01 0.001 0.0001 0.00001
do awk -v var="$i" -F '\t' '($16<var)' chludwig_Y150505_018_IDA_Lib18-Zhonnge_Ecoli_Lib_18.tsv | cut -f11 | sort -u | egrep -i "gms|balrog|smorf|prodigal" | grep -v "XXX" | wc -l
done

cat *tsv | egrep -i "gms|balrog|smorf|prodigal" | grep -v "XXX" > all_canonical_PSMs.tsv
cat *tsv | egrep "XXX" | cut -f10 | cut -f2- -d '.' | rev | cut -f2- -d '.' | rev | sort -u | grep -o -F -f - all_canonical_PSMs.tsv | sort -u | grep -F -f - *tsv | egrep "XXX" | cut -f11 | sort -u > exclude



for j in $(ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/*tsv | rev | cut -f1 -d '/' | rev)
do
for i in 0.01 0.001 0.0001 0.00001 0.000001
do
egrep -i "gms|balrog|smorf|prodigal" /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI*/$j | grep -v "XXX" > temp



egrep "NC_|XXX" /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI*/$j | grep -v -F -f exclude | awk -v var=$i -F '\t' '($16<var)' | cut -f11 | sort -u | grep -v -c "XXX" | sed "s/$/\t$i\tnoncanonical\t$j/g" >> noncan_decoy_counts.tsv
egrep "NC_|XXX" /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI*/$j | grep -v -F -f exclude | awk -v var=$i -F '\t' '($16<var)' | cut -f11 | sort -u | grep -c "XXX" | sed "s/$/\t$i\tdecoy\t$j/g" >> noncan_decoy_counts.tsv
done
done

