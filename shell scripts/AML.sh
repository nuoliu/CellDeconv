#!/bin/sh

#change below to current working directory
cd "/home/nl2/SMART_2019/CellDeconv"

for SAMP_SIZE in 7 9 11 
do
    for alpha in 0.4 0.6 0.8 1.0 1.2 1.4 1.6
    do
        #change the path below
        EXPATH="/home/nl2/SMART_2019/CellDeconv/AML_experiment/sample_${SAMP_SIZE}_alpha_$alpha"
        mkdir $EXPATH
        echo "Sample size=$SAMP_SIZE, alpha=$alpha"
        for i in {1..5}
        do
            
            python ../code/mix_AML.py -o $EXPATH/ -s $SAMP_SIZE
            python lp.py -i $EXPATH/AML_mixture.xlsx -o $EXPATH/ -f True -a $alpha
            python ../code/goodness.py $EXPATH/AML_ClusterMeans.xlsx $EXPATH/AML_mixture_LP_cellTyle_expressions.xlsx $EXPATH/AML_mixture.xlsx $EXPATH/

        done
        rm -r $EXPATH
    done
    
done
