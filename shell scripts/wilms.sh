#!/bin/sh


cd "/home/nl2/SMART_2019/CellDeconv"
EXPATH="/home/nl2/SMART_2019/CellDeconv/data/wilms"
OUT="/home/nl2/SMART_2019/CellDeconv/wilms_experiment"
for tumor in "tumor_118" "tumor_163" "tumor_565"
    do 
    mkdir $OUT/$tumor
        for alpha in 0.4 0.6 0.8 1.0 1.2 1.4 1.6
            do
            current="$OUT/$tumor/alpha_$alpha"
            mkdir current

            # echo "python lp.py -i $EXPATH/sc_mixture.xlsx -o $EXPATH/ -f True -a $alpha"
            echo "$tumor, $alpha"
            python lp_2.py -i $EXPATH/${tumor}.xlsx -o $current -f True -a $alpha
           

            done
    done