#/bin/bash

#author: Annkatrin Goepfert
#Task: Create BED files from Gheerts pipeline for Mariekes mesothelioma CNV files

#Input files from server
##exome
#/home/kopdebeeck/NGS/Mesotheliomen_Marieke/Raw_data_Exomes/CNV_Results
##lpfgrun1
#/home/kopdebeeck/NGS/Mesotheliomen_Marieke/Raw_data_LPFG/CNV_Results
##lpfgrun2
#/home/kopdebeeck/NGS/Mesotheliomen_Marieke/Raw_data_LPFG2/CNV_Results


original_files_exome=/home/agoepfert/documents/mesothelioma_marieke/new_approach/LPWG2/1.calculate_Log_per_bin

#original_filesPath_LPFG=/home/agoepfert/documents/mesothelioma_marieke/original_input/CNV_Results_LPFG2/

output=/home/agoepfert/documents/mesothelioma_marieke/new_approach/LPWG2/3.Rerun_multiIntersect_1-claculate-log-per-bin



cd $output

mkdir -p $output/bed_files

##find input files
find $original_files_exome -type f -name '*.vs.*.txt' | while read line; do
    ##save filename
    filename=`basename $line .txt`
    #echo "Processing file '$line'"
    ##create bed files with 3 columns
    cut -f1,2,3,4 $line > ./$filename.bed
    #cut -f2,3,4,5,6 $line > ./$filename.bed
    sed -i '1d' ./$filename.bed
    sed -n -i '/NaN/!p' ./$filename.bed
    sed -n -i '/Inf/!p' ./$filename.bed
    awk -v OFS="\t\t" '$1=$1' ./$filename.bed > ./$filename.txt #&& cat -vet filename.2.txt
    #split file in two with log is + or -
    #-F defines your field separator, and " " inserts a space. " " -> space, "/t" ->tab
    awk -F" " '{if($4<0)print > "minusLog-'$filename'.bed";else print > "plusLog-'$filename'.bed"}' ./$filename.txt
    mv $filename.bed ./bed_files
    mv $filename.txt ./bed_files
done


#Theshold for plus log files
mkdir -p $output/plusLog

for i in plusLog*.bed
    do
        echo $i
        mv $i ./plusLog
    done



#Theshold for minus log files
mkdir -p $output/minusLog

for i in minusLog*.bed
    do
        echo $i
        mv $i ./minusLog
    done



