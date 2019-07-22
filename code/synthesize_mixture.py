"""
Generate synthetic data for deconvolution
Author:Ivy Liu
"""
import sys
import getopt
import pandas as pd
import numpy as np
import random
def split(cell_exp, depth):
    noise=[(random.random()*2-1)/float(depth) for i in range(len(cell_exp))]
    add_noise=np.array([1.0+x for x in noise])
    subtract_noise=np.array([1.0-x for x in noise])
    return (np.multiply(cell_exp, add_noise), np.multiply(cell_exp,subtract_noise))
def produce_mixture(inputfile=None, output_dir=None, num_samples=None,num_gene=None, num_cells=None):
    cell_names=None
    gene_list=None
    if inputfile:
        data=pd.read_excel(inputfile)
        cell_names=data.columns.values
        if cell_names[0]=="geneSymbol":
            gene_list=data.iloc[:,0].values
            cell_names=cell_names[1:]
    
        num_gene=data.shape[0]
        num_cells=data.shape[1]
        cell_exp=data.T.values   #(#cell, #gene)

    else: #generate synthetic cell line data
        cell_exp=[]
        cell_names=['cell '+str(i+1) for i in range(num_cells)]
        root=np.random.random(num_gene)*400     #scaled doesn't really matter
        cell_exp.append(root)
        num_cells_made=1
        depth=1
        while True:
            if num_cells==num_cells_made:
                break
            new_cells=[]
            if num_cells>=2*num_cells_made:  #need to split everyone we have and more...
                for old_cell in cell_exp:
                    new1,new2=split(old_cell,depth)
       
                    new_cells.append(new1)
                    new_cells.append(new2)
                cell_exp=new_cells
                num_cells_made=len(cell_exp)
                depth+=1
            else:  #num_cells>num_cells_made but less than twice
                num_more=num_cells-num_cells_made
                #split the front 
                for i in range(num_more):
                    new1,new2=split(cell_exp[i],depth)
       
                    new_cells.append(new1)
                    new_cells.append(new2)
                #keep the old ones on the back
                for i in range(num_more,num_cells_made):
  
                    new_cells.append(cell_exp[i])
                cell_exp=new_cells
                break
        
        gene_list=["Gene "+str(i+1) for i in range(num_gene)]
        
        #write out the generated cell lines
        out_cellLines=np.array(cell_exp)
        out_cellLines=out_cellLines.T
        out_cellLines=pd.DataFrame(out_cellLines,columns=cell_names)
        out_cellLines.insert(0,"geneSymbol",gene_list)
        cell_fn='synthetic_cellLines.xlsx'
        if output_dir:
            cell_fn=output_dir+cell_fn
        with pd.ExcelWriter(cell_fn) as writer:
            out_cellLines.to_excel(writer, index=False)



    composition=np.zeros((num_samples,num_cells))
        #randomly generate the compositions:
    for i in range(num_samples):
        dist=np.random.random(num_cells)
        dist/=dist.sum()
        composition[i]=dist
    

    mixture=np.zeros((num_samples,num_gene))
    for i in range(num_samples):
        for j in range(num_cells):
            mixture[i]+=composition[i][j]*cell_exp[j]

    mixture=mixture.T
    out_mixture=pd.DataFrame(mixture,columns=["sample "+str(i+1) for i in range(num_samples)])
    out_mixture.insert(0,'geneSymbol',gene_list)
    mixture_fn="synthetic_mixture.xlsx"


    composition=composition.T  #(#celss, #sample)
    out_compos=pd.DataFrame(composition, columns=["sample "+str(i+1) for i in range(num_samples)])
    out_compos.insert(0,'cell line',cell_names)
    compo_fn="synthetic_composition.xlsx"
    if output_dir:
        compo_fn=output_dir+compo_fn
        mixture_fn=output_dir+mixture_fn

    with pd.ExcelWriter(compo_fn) as writer:
        out_compos.to_excel(writer, index=False)

    with pd.ExcelWriter(mixture_fn) as writer:
        out_mixture.to_excel(writer, index=False)
def main(argv):
    inputfile = None
    outputfile = None
    num_sample=0
    num_gene=0
    num_cells=0
    try:
        opts, args = getopt.getopt(argv,"hi:o:s:")
    except getopt.GetoptError:
        print('synthesize_mixture.py -i <input excel file> -o <out_directory> -s num_samples or synthesize_mixture.py -i num_gene,num_cell -o <out_directory> -s num_samples')
        sys.exit(2)
    for opt, arg in opts:
       
        if opt == '-h':
            print('synthesize_mixture.py -i <input excel file> -o <out_directory> -s num_samples or synthesize_mixture.py -i num_gene,num_cell -o <out_directory> -s num_samples')
            sys.exit()
        elif opt=='-i':
            info=arg.split('.')
            if len(info)>1:
                inputfile = arg
            else:
                info=arg.split(',')
                num_gene=int(info[0])
                num_cells=int(info[1])
        elif opt=="-o":
            outputfile = arg
        elif opt=="-s":
            num_sample=int(arg)
        


    if inputfile!=None:
        # print("Generate synthetic mixture data using cell lines from %s"%inputfile)
        produce_mixture(inputfile,outputfile,num_sample, None, None)
    else:
        # print("Generate synthetic mixture data randomly...")
        produce_mixture(None, outputfile, num_sample,num_gene, num_cells)




if __name__ == "__main__":
   main(sys.argv[1:])