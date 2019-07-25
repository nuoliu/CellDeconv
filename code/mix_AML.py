import sys
import getopt
import pandas as pd
import numpy as np
import random
from sklearn.preprocessing import normalize
from scipy.stats import variation
def produce_mixture(out_dir, num_sample):

    expression_data=pd.read_excel("../data/AML/AML328-D113_expression.xlsx")
    clustered_sample=pd.read_excel("../data/AML/AML328-D113_selectedSamples.xlsx")
    num_genes=expression_data.values.shape[0]
    gene_list=expression_data.loc[:,"Gene"].values
    clusterMembers={}
    clusters=clustered_sample["CellType"].unique()
    num_clusters=len(clusters)
    cluster_sizes=np.zeros(num_clusters)
    for index,cluster in enumerate(clusters):
        cluster_member=clustered_sample.loc[clustered_sample["CellType"]==cluster]
        cluster_member=cluster_member["Cell"].values.tolist()
        cluster_sizes[index]=len(cluster_member)
        clusterMembers[cluster]=cluster_member

    #start with composistino of raw counts
    #indicates how many barcodes to select from each cluster
    compositions=np.zeros((num_sample,num_clusters),dtype=np.int)

    expressions=np.zeros((num_sample,num_genes))
    for i in range(num_sample):
        ingredients=[]
        for j, cluster in enumerate(clusters):
            compositions[i][j]=random.randint(1,cluster_sizes[j])
            #for each sample, fetch the barcodes indices to be mixed 
            sample_from_cluster=random.sample(clusterMembers[cluster],compositions[i][j])
            ingredients.extend([expression_data.loc[:,x].values for x in sample_from_cluster])   #add this list of arrays
        #mix/average all the ingredients (expression from each barcodes selected)
        ingredients=np.array(ingredients)
        sample_expression=np.mean(ingredients,axis=0)
        expressions[i]=sample_expression

    #convert raw composition into ratios
    normed_composition = normalize(compositions, axis=1, norm='l1')

    #calculate ground truth
    cluster_means=np.zeros((num_clusters,num_genes))
    for j, cluster in enumerate(clusters):
        all_members=clusterMembers[cluster]
        all_member_expression=np.array([expression_data.loc[:,x].values for x in all_members])
        cluster_mean=np.mean(all_member_expression,axis=0)
        cluster_means[j]=cluster_mean
   
    





    #prepare output dataframs
    normed_composition=normed_composition.T  #(#cluster, #sample)
    out_compos=pd.DataFrame(normed_composition, columns=["sample "+str(i+1) for i in range(num_sample)])
    out_compos.insert(0,'cluster',["cluster"+str(j) for j in range(num_clusters)])
    compo_fn=out_dir+"AML_composition.xlsx"

    expressions=expressions.T
    out_mixture=pd.DataFrame(expressions, columns=["sample "+str(i+1) for i in range(num_sample)])
    out_mixture.insert(0,'geneSymbol',gene_list)
    #filter genes
    #filter out genes with max lower than 1
    MAX=np.amax(expressions, axis=1)
    out_mixture=out_mixture.loc[MAX>0.5]
    expressions=out_mixture.iloc[:,1:].values
    CV=variation(expressions, axis=1)
    out_mixture=out_mixture.loc[CV>0.5]
    out_gene=out_mixture["geneSymbol"].values

    mixture_fn=out_dir+"AML_mixture.xlsx"

    
    cluster_means=cluster_means.T
    out_clusterMeans=pd.DataFrame(cluster_means,columns=["cluster"+str(j) for j in range(num_clusters)])
    out_clusterMeans.insert(0,"geneSymbol",gene_list)
    out_clusterMeans=out_clusterMeans.loc[out_clusterMeans["geneSymbol"].isin(out_gene)]
    clusterMean_fn=out_dir+"AML_ClusterMeans.xlsx"
    #output files: composistion, ground truth (cluster means), mixture
    with pd.ExcelWriter(compo_fn) as writer:
        out_compos.to_excel(writer, index=False)
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
        print('mix_AML.py -o <out_directory> -s num_samples')
        sys.exit(2)
    for opt, arg in opts:
  
        if opt == '-h':
            print('mix_AML.py -o <out_directory> -s num_samples')
            sys.exit()
        elif opt=="-o":
            out_dir = arg
        elif opt=="-s":

            num_sample=int(arg)
        


    
    # print("Generate synthetic mixture from single cell RNA data...")
    produce_mixture(out_dir, num_sample)




if __name__ == "__main__":
   main(sys.argv[1:])