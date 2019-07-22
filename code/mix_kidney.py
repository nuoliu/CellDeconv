
import pandas as pd
import random
import sys
import getopt
import numpy as np
from scipy.stats import variation
def produce_mixture(input_file, num_sample, num_cells, out_dir):

    data=pd.read_excel(input_file)
    selected_cells=random.sample(list(range(1,data.shape[1])),num_cells)
    selected_cell_names=[data.columns.values[x] for x in selected_cells]
    #select genes
    expressions=data.iloc[:,1:].values
    CV=variation(expressions, axis=1)
    data=data.loc[CV>1.0]
    gene_list=data.iloc[:,0].values
    num_gene=len(gene_list)
    #write out the true cell lines
    out_l=[0]
    out_l.extend(selected_cells)
    trueCell=data.iloc[:,out_l]
    cell_fn=out_dir+'kidney_cellLines.xlsx'
    with pd.ExcelWriter(cell_fn) as writer:
        trueCell.to_excel(writer, index=False)

    composition=np.zeros((num_sample,num_cells))
        #randomly generate the compositions:
    l=[0.5,0.25, 0.15, 0.05, 0.05]
    for i in range(num_sample):
        # dist=np.random.random(num_cells)
        # dist/=dist.sum()
        # composition[i]=dist
        random.shuffle(l)
        composition[i]=np.array(l)

    
  

    mixture=np.zeros((num_sample,num_gene))
  
    for i in range(num_sample):
        for j in range(num_cells):
            mixture[i]+=composition[i][j]*data.iloc[:,selected_cells[j]].values

    mixture=mixture.T
    out_mixture=pd.DataFrame(mixture,columns=["sample "+str(i+1) for i in range(num_sample)])
    out_mixture.insert(0,'geneSymbol',gene_list)
    mixture_fn=out_dir+"kidney_mixture.xlsx"

    composition=composition.T  #(#celss, #sample)
    out_compos=pd.DataFrame(composition, columns=["sample "+str(i+1) for i in range(num_sample)])
    out_compos.insert(0,'cell line',selected_cell_names)
    compo_fn=out_dir+"kidney_composition.xlsx"

    with pd.ExcelWriter(compo_fn) as writer:
        out_compos.to_excel(writer, index=False)

    with pd.ExcelWriter(mixture_fn) as writer:
        out_mixture.to_excel(writer, index=False)


def main(argv):  
    
    out_dir = None
    inputfile=sys.argv[1]
    num_sample=int(sys.argv[2])
    num_cells=int(sys.argv[3])
    out_dir=sys.argv[4]


    
    # print("Generate synthetic mixture from single cell RNA data...")
    produce_mixture(inputfile,num_sample, num_cells,out_dir)




if __name__ == "__main__":
   main(sys.argv[1:])