#Mori2021

cut -f2 promising | sort -u | cut -f1 -d "(" | grep -v "NC_000913.3_101381" | grep -v "NC_000913.3_4041" | while IFS= read -r i; do
    grep -F "$i(" promising | cut -f1 -d ":" | rev | cut -f1 -d '/' | cut -f2- -d "." | rev | while IFS= read -r j; do
        grep -F "$i(" "/stor/scratch/Ochman/hassan/112724_protogene_extension/data/Mori2021/MS/mgf/real_mgf/$j.tsv" | 
        awk -F '\t' '($16<0.00001)' | cut -f4 | rev | cut -f2 -d "\"" | rev | sed 's/^/"/g' | sed 's/$/"/g' | while IFS= read -r k; do
            sed -n "/$k/,/END IONS/p" "${j}_annotated.mgf" > "${i}_${j}.spectra"
        done
    done
done
