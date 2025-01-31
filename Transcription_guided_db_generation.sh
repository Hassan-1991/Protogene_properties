#Caglar2017:

awk -F '\t' '{OFS=""}{print "time htseq-count -f bam -a 0 -t CDS --nonunique all /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/RNAseq\/bamfiles\/",$1,"_sorted_pairedonly_byname.bam REL606.final.gtf \| head -n -5 \| sed \"s\/^\/",$1,"\t\/g\" >> ",$1,"_htseq.tsv"}' /stor/work/Ochman/hassan/proteomics_denovo/1031_RNAseq_databasemaking/MURIfiles_with_RNAprotein_data.txt | sed "s/\.mgf//g" > running.sh

awk -F '\t' '{print $5-$4,$9}' REL606.final.gtf | cut -f 1 -d ';' | sed "s/transcript id \"//g" | sed "s/\"//g" | sed "s/ transcript_id /\t/g" | awk -F '\t' '{OFS=FS}{print $2,$1}' | sort -k1 > ORF_lengths.txt

for i in MURI*htseq.tsv
do
sort -k2 $i | join -1 2 -2 1 - ORF_lengths.txt | sed "s/ /\t/g" | awk -F '\t' '{OFS=FS}{print $1, $2, $3, $4, $3/$4}' > interim
sum=$(awk '{ sum += $5 } END { print sum }' interim)
awk -v sum="$sum" 'BEGIN{OFS="\t"} {print $0, ($5/sum) * 1000000}' interim > $i.tpm
awk '($6>0.5)' $i.tpm | cut -f 1 | sort -u | sed "s/.*/\"&\"/g" > $i.90quart
grep -F -f $i.90quart REL606.final.gtf | gtf2bed | bedtools getfasta -s -name -fi ../REL606.faa -bed - > $i.database
/stor/work/Ochman/hassan/tools/faTrans -stop $i.database $i.database.proteins.fa
done

#Total proteins: 73265
#Canonical: 3474
#Non-canonical: 69791
#Median size of database:
for i in *tsv.tpm; do awk -F '\t' '($6>0.5)' $i | cut -f1 | sort -u | wc -l; done | sort -n | awk '{data[NR] = $1} END {if (NR % 2 == 1) print data[int(NR/2) + 1]; else print (data[int(NR/2)] + data[int(NR/2) + 1]) / 2}'

#generate mzid files:
cat /stor/work/Ochman/hassan/proteomics_denovo/1031_RNAseq_databasemaking/MURIfiles_with_RNAprotein_data.txt | awk '{OFS=""}{print "time java -Xmx3500M -jar \/stor\/work\/Ochman\/hassan\/Fall_2022\/massspec_database\/MSGFplus\/MSGFPlus.jar -s /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/",$1," -d /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/",$1,"_htseq.tsv.database.proteins.fa -inst 1 -t 10ppm -ti 0,1 -mod \/stor\/work\/Ochman\/hassan\/Fall_2022\/08112023_MSGFpercolator\/mods.txt -ntt 2 -tda 1 -maxMissedCleavages 2 -addFeatures 1 -thread 72"}' > running.sh

#convert to tsv:
cd /stor/work/Ochman/hassan/Fall_2022/massspec_database/MSGFplus/ #For some reason this code only runs in the directory of the jar file
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/*mzid | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i "$1" -o "$1".tsv -showQValue 1 -showDecoy 1 -unroll 1"}' | sed "s/mzid.tsv/tsv/g" > running.sh

for i in $(ls /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/*mzid | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
do
echo "python3 /stor/work/Ochman/hassan/mass_spec/test_files/mgf_search_result_annotator_test1.py --format MSGF_ident --input /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/MURI*/${i}.mgf --search /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/${i}.mzid --output /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/${i}_annotated.mgf"
done > running.sh

cat MURI*tsv | grep "REL606" | grep -v "XXX" > all_canonical_PSMs.tsv
cat MURI*tsv | egrep "XXX|NC_" | cut -f10 | cut -f2- -d '.' | rev | cut -f2- -d '.' | rev | sort -u | grep -o -F -f - all_canonical_PSMs.tsv | sort -u | grep -F -f - *tsv | egrep "XXX|NC_" | cut -f11 | sort -u > exclude
grep -v -F -f exclude MURI*tsv | awk -F '\t' '($16<0.00001)' | cut -f11 | sort -u > test2
grep -v "REL606" test2 | grep -F -f - MURI*tsv | awk -F '\t' '($16<0.00001)' | cut -f1,11 | grep -v "XXX" > promising_spectra

cut -f2 promising_spectra | sort -u | cut -f1 -d "(" | grep -v "NC_012967.1_102979" | grep -v "NC_012967.1_122979" | while IFS= read -r i; do
    grep -F "$i(" promising_spectra | cut -f1 -d ":" | rev | cut -f2- -d "." | rev | while IFS= read -r j; do
        grep -F "$i(" "/stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/$j.tsv" | 
        awk -F '\t' '($16<0.00001)' | cut -f4 | rev | cut -f2 -d "\"" | rev | while IFS= read -r k; do
            sed -n "/$k$/,/END IONS/p" "${j}_annotated.mgf" > "${i}_${j}.spectra"
        done
    done
done

