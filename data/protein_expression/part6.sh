#!/bin/sh


cd "/home/nl2/SMART_2019/CellDeconv"
EXPATH="/home/nl2/SMART_2019/CellDeconv/protein_expression"
for patient in "patient6"  
    do 
    filename="${patient}.xlsx"
    mkdir $EXPATH/$patient
        for alpha in 0.4 0.6 0.8 1.0 1.2 1.4 1.6
            do
            current="$EXPATH/$patient/$alpha/"
            mkdir current

            # echo "python lp.py -i $EXPATH/sc_mixture.xlsx -o $EXPATH/ -f True -a $alpha"
            echo "$patient, $alpha"
            python lp_2.py -i $EXPATH/$filename -o $current -f True -a $alpha
           

            done
    done