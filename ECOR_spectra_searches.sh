#MSGF Similar to running other searches
#Starting from mzid:
cd /stor/work/Ochman/hassan/Fall_2022/massspec_database/MSGFplus/ #For some reason this code only runs in the directory of the jar file
ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/*mzid | awk '{OFS=""}{print "java -Xmx3500M -cp MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i "$1" -o "$1".tsv -showQValue 1 -showDecoy 1 -unroll 1"}' | sed "s/mzid.tsv/tsv/g" > running.sh

for i in $(ls /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/*mzid | rev | cut -f 1 -d '/' | rev | cut -f 1 -d '.')
do
echo "python3 /stor/work/Ochman/hassan/mass_spec/test_files/mgf_search_result_annotator_test1.py --format MSGF_ident --input /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/${i}.mgf --search /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/${i}.mzid --output /stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/${i}_annotated.mgf"
done > running.sh

cat 8925*tsv | egrep -v "XXX" | egrep -i "prodigal|gms|balrog|smorf" > all_canonical_PSMs.tsv
cat 8925*tsv | egrep "XXX" | cut -f10 | cut -f2- -d '.' | rev | cut -f2- -d '.' | rev | sort -u | grep -o -F -f - all_canonical_PSMs.tsv | sort -u | grep -F -f - *tsv | egrep "XXX" | cut -f11 | sort -u > exclude
cat 8925*tsv | egrep -iv "XXX|gms|balrog|smorf|prodigal" | cut -f10 | cut -f2- -d '.' | rev | cut -f2- -d '.' | rev | sort -u | grep -o -F -f - all_canonical_PSMs.tsv | sort -u | grep -F -f - *tsv | egrep -iv "XXX|gms|balrog|smorf|prodigal" | cut -f11 | sort -u >> exclude
grep -v -F -f exclude 8925*tsv | awk -F '\t' '($16<0.0001)' | cut -f11 | sort -u > test2
grep -v "REL606" test2 | grep -F -f - MURI*tsv | awk -F '\t' '($16<0.00001)' | cut -f1,11 | grep -v "XXX" > promising_spectra

egrep -vi "gms|prodigal|smorf|balrog|XXX" test2 | grep -F -f - 8925*tsv | awk -F '\t' '($16<0.00001)' | cut -f1,11 | grep -v "XXX" > promising_spectra

cut -f2 promising_spectra | sort -u | cut -f1 -d "(" | while IFS= read -r i; do
    grep -F "$i(" promising_spectra | cut -f1 -d ":" | rev | cut -f2- -d "." | rev | while IFS= read -r j; do
        grep -F "$i(" "/stor/scratch/Ochman/hassan/112724_protogene_extension/data/ECOR_2023/$j.tsv" | 
        awk -F '\t' '($16<0.00001)' | cut -f4 | rev | cut -f2 -d "\"" | rev | while IFS= read -r k; do
            sed -n "/$k$/,/END IONS/p" "${j}_annotated.mgf" > "${i}_${j}.spectra"
        done
    done
done

