import sys
import getopt
import pandas as pd
import numpy as np
import random
from sklearn.preprocessing import normalize

def readClusters():
    fn=open("data/scRNA/Filtered_clusters.txt")
    clusters={}
    for line in fn:
        info=line.rstrip().split()
        clusters[info[0][-1]]=[int(x) for x in info[1:]]

    return clusters
def produce_mixture(out_dir, num_sample=1):

    expression_data=pd.read_excel("data/scRNA/formatted_expression.xlsx")
    num_genes=expression_data.values.shape[0]
    gene_list=expression_data.loc[:,"geneSymbol"].values
    clusterMembers=readClusters()
    num_clusters=len(clusterMembers)
    cluster_sizes=np.zeros(num_clusters)
    for key, value in clusterMembers.items():
        cluster_sizes[int(key)]=len(value)
    
    num_samples=expression_data.shape[1]-1
    # expressions=expression_data.loc[:,random.randint(0,num_samples)].values
    j=random.randint(0,num_clusters-1)
    sample_from_cluster=random.sample(clusterMembers[str(j)],random.randint(1,cluster_sizes[j]))
    ingredients=[]
    ingredients.extend([expression_data.loc[:,x].values for x in sample_from_cluster]) 
    ingredients=np.array(ingredients)
    expressions=np.mean(ingredients,axis=0)
   


    #calculate ground truth
    cluster_means=np.zeros((num_clusters,num_genes))
    for j in range(num_clusters):
        all_members=clusterMembers[str(j)]
        all_member_expression=np.array([expression_data.loc[:,x].values for x in all_members])
        cluster_mean=np.mean(all_member_expression,axis=0)
        cluster_means[j]=cluster_mean
    


    expressions=expressions.T
    out_mixture=pd.DataFrame(expressions, columns=["sample 1"])
    out_mixture.insert(0,'geneSymbol',gene_list)
    mixture_fn=out_dir+"sc_mixture.xlsx"

    cluster_means=cluster_means.T
    out_clusterMeans=pd.DataFrame(cluster_means,columns=["cluster"+str(j) for j in range(num_clusters)])
    out_clusterMeans.insert(0,"geneSymbol",gene_list)
    clusterMean_fn=out_dir+"sc_ClusterMeans.xlsx"


    #output files: composistion, ground truth (cluster means), mixture
 
    with pd.ExcelWriter(mixture_fn) as writer:
        out_mixture.to_excel(writer, index=False)
    with pd.ExcelWriter(clusterMean_fn) as writer:
        out_clusterMeans.to_excel(writer, index=False)
    

def main(argv):  
    
    out_dir = None
    num_sample=0

    try:
        opts, args = getopt.getopt(argv,"hs:o:")
    except getopt.GetoptError:
        print('mix_sc.py -o <out_directory> -s num_samples')
        sys.exit(2)
    for opt, arg in opts:
  
        if opt == '-h':
            print('mix_sc.py -o <out_directory> -s num_samples')
            sys.exit()
        elif opt=="-o":
            out_dir = arg
        elif opt=="-s":

            num_sample=int(arg)
        


    
    # print("Generate synthetic mixture from single cell RNA data...")
    produce_mixture(out_dir, num_sample)




if __name__ == "__main__":
   main(sys.argv[1:])