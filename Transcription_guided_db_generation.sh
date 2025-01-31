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

#custom MSGF code:
cat /stor/work/Ochman/hassan/proteomics_denovo/1031_RNAseq_databasemaking/MURIfiles_with_RNAprotein_data.txt | awk '{OFS=""}{print "time java -Xmx3500M -jar \/stor\/work\/Ochman\/hassan\/Fall_2022\/massspec_database\/MSGFplus\/MSGFPlus.jar -s /stor/scratch/Ochman/hassan/112724_protogene_extension/data/Caglar2017/MS_mgf/",$1," -d /stor/scratch/Ochman/hassan/112724_protogene_extension/RNAseq_diversion/",$1,"_htseq.tsv.database.proteins.fa -inst 1 -t 10ppm -ti 0,1 -mod \/stor\/work\/Ochman\/hassan\/Fall_2022\/08112023_MSGFpercolator\/mods.txt -ntt 2 -tda 1 -maxMissedCleavages 2 -addFeatures 1 -thread 72"}' > running.sh
