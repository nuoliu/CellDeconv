"""
interprets the result of deconvolution and 
measures the goodness of prediction
Author: Ivy Liu
Summer 2019
"""
import sys
import pandas as pd
import numpy as np
from math import exp
from scipy.stats.stats import pearsonr
import networkx as nx
from synthesize_mixture import produce_mixture

def calculateMetric(truth_fn,inferred_fn, mixture_fn):
    def convertRtoZ(r):
        if r==1:
            return 10000
        else:
            return 0.5*np.log((1.0+r)/float(1.0-r))

    true_cells=pd.read_excel(truth_fn)
    num_gene=true_cells.values.shape[0]
    true_cells=true_cells.drop(columns="geneSymbol")
    true_data=true_cells.values.T
    inferred_cells=pd.read_excel(inferred_fn)
    inferred_cells=inferred_cells.drop(columns="geneSymbol")
    infer_data=inferred_cells.values.T
    assert (true_data.shape[1]==infer_data.shape[1]),"Cell line dimension don't match! "
    num_tru=true_data.shape[0]
    num_inferred=infer_data.shape[0]
    mixture_samples=pd.read_excel(mixture_fn).drop(columns="geneSymbol")
    mixture_data=mixture_samples.values.T
    dummpy_exp=np.mean(mixture_data, axis = 0)   #average of all true cell lines
  
    G=nx.Graph()
    #TODO calculate correlation coefficients and assign greedily
    dimension=max(num_inferred,num_tru)
    corr=np.zeros((dimension,dimension))
    
    #use dummy_exp for the mismatch in dimension
    #add edges to the undirected graph
    for i in range(dimension):
        for j in range(dimension):
            if i>=num_tru or j>=num_inferred:
                #add dummy nodes in bipartite graph 
                G.add_edge("inferred "+str(j),"True "+str(i),weight=0.0)
               
            else:
                G.add_edge("inferred "+str(j),"True "+str(i),weight=pearsonr(infer_data[j],true_data[i])[0])


   #use maximum weighted matching algorithm to find 
   # the maximal matching that maximizes the sum of weights(correlation coefficient)
    max_match=nx.max_weight_matching(G,maxcardinality=True) #This function takes time O(number_of_nodes ** 3)

    z_scores=[]
    for u, v in max_match:
        if "inferred" in u:
            inferred_idx=int(u.split()[1])
            true_idx=int(v.split()[1])
        else:
            inferred_idx=int(v.split()[1])
            true_idx=int(u.split()[1])
        #if it was originally matched to dummy vertex with zero weight on edge
        #we penalize the extra vertex on either side by allowing it
        #to match to the dummy expression which is the average of true cell lines
        if inferred_idx>=num_inferred:
            r=pearsonr(dummpy_exp,true_data[true_idx])[0]
        elif  true_idx>=num_tru:
            r=pearsonr(dummpy_exp,infer_data[inferred_idx])[0]
        else:
            r=G[u][v]["weight"]

        z_scores.append(convertRtoZ(r))
       
    z_avg=sum(z_scores)/len(z_scores)
    r_avg=(exp(2*z_avg)-1.0)/(exp(2*z_avg)+1.0)

    #calculate what would be outputted if no cell line is inferred
    z_dum=0
    for i in range(num_tru):
        z_dum+=convertRtoZ(pearsonr(dummpy_exp,true_data[i])[0])
    z_dum/=num_tru
    r_dum=(exp(2*z_dum)-1.0)/(exp(2*z_dum)+1.0)
    
    return (num_inferred, num_tru, r_avg, r_dum, num_gene)


def main():
    if len(sys.argv)<5:
        print("wrong number of arguments, please provide both true cell line and inferred cell lines!")
        sys.exit(1)
    else:
        truth_fn=sys.argv[1]
        inferred_fn=sys.argv[2]
        mixture_fn=sys.argv[3]
        out_dir=sys.argv[4]
        num_inferred, num_tru,r_infer, r_dum, num_gene=calculateMetric(truth_fn,inferred_fn,mixture_fn)
        # print("Matching %d inferred cell lines to %d true cell lines..."%(num_inferred,num_tru))
        # print("r value of inferrence is %5.4f"%r_avg)
        # print("r value of dummy inferrence is %5.4f"%r_dum)
        r_value=0
        for i in range(5):
            produce_mixture(None, out_dir, 1,num_gene, num_inferred)
            random_cellLine=out_dir+'synthetic_cellLines.xlsx'
            _, _, r_avg,_, _=calculateMetric(truth_fn,random_cellLine,mixture_fn)
            r_value+=r_avg
        r_value/=5
        # print("The average r value from %d randomly synthesized cell lines is %5.4f" %(num_inferred,r_value))
        print("%d,%5.4f,%5.4f,%5.4f"%(num_inferred,r_infer,r_dum,r_value))


   
       



if __name__ == "__main__":
    main()