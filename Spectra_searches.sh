cd /stor/scratch/Ochman/hassan/112724_protogene_extension

awk '{OFS=""}{print "time java -Xmx3500M -jar \/stor\/work\/Ochman\/hassan\/Fall_2022\/massspec_database\/MSGFplus\/MSGFPlus.jar -s \/stor\/work\/Ochman\/hassan\/mass_spec\/raw_data\/",$1," -d \/stor\/work\/Ochman\/hassan\/proteomics_denovo\/1031_RNAseq_databasemaking\/",$1,"_htseq.tsv.database.tpm1.proteins.fa -inst 1 -t 10ppm -ti 0,1 -mod \/stor\/work\/Ochman\/hassan\/Fall_2022\/08112023_MSGFpercolator\/mods.txt -ntt 2 -tda 1 -maxMissedCleavages 2 -addFeatures 1 -thread 72"}' MURIfiles_with_RNAprotein_data.txt | sed "s/\.mgf//g" > MSGF.sh

#Convert to tsv
cd /stor/work/Ochman/hassan/Fall_2022/massspec_database/MSGFplus/ #For some reason this code only runs in the directory of the jar file
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI_*/*mzid | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i "$1" -o "$1".tsv -showQValue 1 -showDecoy 1 -unroll 1"}' | sed "s/mzid.tsv/tsv/g" > mzid2tsv.sh


ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI_*/*mzid | rev | cut -f 1 -d '/' | cut -f 2- -d '.' | rev | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i /stor/work/Ochman/hassan/proteomics_denovo/htseq/renew/",$1," -o /stor/work/Ochman/hassan/proteomics_denovo/htseq/renew/",$1,".tsv -showQValue 1 -showDecoy 1 -unroll 1"}' > mzidtotsv.sh

ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI_*/*mzid | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i /stor/work/Ochman/hassan/proteomics_denovo/htseq/renew/",$1," -o /stor/work/Ochman/hassan/proteomics_denovo/htseq/renew/",$1,".tsv -showQValue 1 -showDecoy 1 -unroll 1"}'
#ms2rescore
