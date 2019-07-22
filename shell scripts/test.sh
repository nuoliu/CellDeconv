#!/bin/sh


cd "/home/nl2/SMART_2019/CellDeconv"
EXPATH="/home/nl2/SMART_2019/CellDeconv/experiment"
for i in {1..20}
    do
       
        python mix_kidney_random.py data/kidneyCCL_highCVgenes.xlsx 6 5 $EXPATH/
        python goodness_zero.py $EXPATH/kidney_cellLines.xlsx $EXPATH/kidney_mixture.xlsx

    done