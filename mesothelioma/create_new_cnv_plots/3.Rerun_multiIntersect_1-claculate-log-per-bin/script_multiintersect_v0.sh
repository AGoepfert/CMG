#/bin/bash

#bedfile folder
bedFile=/home/agoepfert/documents/mesothelioma_marieke/new_approach/LPWG2/3.Rerun_multiIntersect_1-claculate-log-per-bin

sampleALL=""
filenameALL=""


for i in $bedFile/minusLog/minusLog*.bed
    do
        echo `basename $i`
        #filename=`$i`
        sample=`basename $i .bed`
        #echo $sample #"_________" $filename
        sampleALL=`printf "%s $sample" "$sampleALL"`
        filenameALL=`printf "%s $i" "$filenameALL"`
        #echo $filenameALL
        #echo $sampleALL
    done
 

multiIntersectBed -i $filenameALL -names $sampleALL > out_minus.txt

sampleALLplus=""
filenameALLplus=""


for j in $bedFile/plusLog/plusLog*.bed
    do
        echo `basename $j`
        #filename=`$j`
        sample=`basename $j .bed`
        #echo $sample #"_________" $filename
        sampleALLplus=`printf "%s $sample" "$sampleALLplus"`
        filenameALLplus=`printf "%s $j" "$filenameALLplus"`
        #echo $filenameALLplus
    done

multiIntersectBed -i $filenameALLplus -names $sampleALLplus > out_plus.txt
