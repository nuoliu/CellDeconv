#!/bin/sh


cd "/home/nl2/SMART_2019/CellDeconv"

for SAMP_SIZE in 11 13 15
do
    for alpha in 0.4 0.6 0.8 1.0 1.2 1.4 1.6
    do
        EXPATH="/home/nl2/SMART_2019/CellDeconv/kidney_experiment_random/sample_${SAMP_SIZE}_alpha_$alpha"
        mkdir $EXPATH
        echo "Sample size=$SAMP_SIZE, alpha=$alpha"
        for i in {1..5}
        do
            
            # echo "python mix_sc.py -o $EXPATH/ -s $SAMP_SZIE" 
            python mix_kidney_random.py data/kidneyCCL_highCVgenes.xlsx $SAMP_SIZE 5 $EXPATH/
            # echo "python lp.py -i $EXPATH/sc_mixture.xlsx -o $EXPATH/ -f True -a $alpha"
            python lp.py -i $EXPATH/kidney_mixture.xlsx -o $EXPATH/ -f True -a $alpha
            # echo "python goodness.py $EXPATH/sc_ClusterMeans.xlsx $EXPATH/sc_mixture_LP_cellTyle_expressions.xlsx $EXPATH/sc_mixture.xlsx $EXPATH/"
            python goodness.py $EXPATH/kidney_cellLines.xlsx $EXPATH/kidney_mixture_LP_cellTyle_expressions.xlsx $EXPATH/kidney_mixture.xlsx $EXPATH/

        done
        rm -r $EXPATH
    done
    
done