import pandas as pd
import numpy as np
from scipy.stats import variation
def formatMatrix():
    # with open("matrix.mtx") as file:
    #     gene=[]
    #     barcode=[]
    #     expression=[]
    #     for line in file:
    #         try:
    #             info=line.split()
    #             gene.append(int(info[0]))
    #             barcode.append(int(info[1]))
    #             expression.append(int(info[2]))
    #         except:
    #             pass
    # data={'geneSymbol':gene, 'barcode':barcode,'expression':expression}
    # data=pd.DataFrame(data)
    data=pd.read_excel("matrix.xlsx")
    num_barcodes=data["barcode"].unique().shape[0]
    #filter out only genes that appear in more than half of the barcodes
    # groups=data.groupby("geneSymbol")
    # gene_counts=data.groupby("geneSymbol").count()
    # genes=gene_counts.loc[gene_counts["barcode"]>=(num_barcodes/2)].index
    # data=data.loc[data["geneSymbol"].isin(genes)]
    print(len(data["barcode"].unique()))

    data=data.groupby("geneSymbol").filter(lambda x: len(x)>=0.3*num_barcodes)
    genes=data["geneSymbol"].unique()
    #get all the barcodes left
    barcodes=data["barcode"].unique()
    print(len(data["barcode"].unique()))

    matrix=np.zeros((genes.shape[0],barcodes.shape[0]))
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            try:
                matrix[i][j]=data.loc[(data['geneSymbol']==genes[i]) & (data['barcode'] == barcodes[j])]['expression']
            except:
                matrix[i][j]=0
    



    # genes=data.groupby('geneSymbol').count().nlargest(num_genes, columns='barcode').index
    # data=data.loc[data["geneSymbol"].isin(genes)]
  


    out_matrix=pd.DataFrame(matrix, columns=barcodes)
    out_matrix.insert(0,"geneSymbol",genes)
    
    

    clusters=pd.read_excel("ClusterMembers.xlsx")
    valid_barcodes=barcodes
    codes=convertBarcodes()
    cluster_names=clusters.columns.values
    barcodes=clusters.values.T #each row is a list of barcodes that belong to the cluster
    barcodes=[ [ x for x in barcode if isinstance(x, str)] for barcode in barcodes] #remove nan
    filtered_clusters=[]
    for cluster in barcodes:
        filtered_cluster=[]
        for barcode in cluster:
            #if the barcode index is in the filtered 342 barcodes, keep it
            numerical=codes[barcode]
            if numerical in valid_barcodes:
                filtered_cluster.append(numerical)
        filtered_clusters.append(filtered_cluster)
    f=open("Filtered_clusters.txt",'w')
    for index, cluster in enumerate(filtered_clusters):
        barcodes=' '.join(str(v) for v in cluster)
        line="cluster"+str(index)+" "+barcodes+"\n"
        f.write(line)
    f.close()
    #in the formatted expression matrix filter out only genes with high variation across clusters
    CVs=[getCV(out_matrix.iloc[x],filtered_clusters) for x in range(out_matrix.shape[0])]
    out_matrix["CV"]=CVs
    out_matrix=out_matrix.loc[out_matrix["CV"]>1.2]
    out_matrix=out_matrix.drop(columns="CV")
    print(out_matrix.shape)
            
    with pd.ExcelWriter("formatted_expression.xlsx") as writer:
        out_matrix.to_excel(writer, index=False)
            


def convertBarcodes():
    f=open('barcodes.tsv')
    codes={}
    index=1
    for line in f:
        code=line.split('-')[0]
        codes[code]=index
        index+=1
    return codes

def getCV(row_data, filtered_clusters):
    cluster_exp=[]
    for cluster in filtered_clusters:
        cluster_exp.append(row_data[cluster].values)

    cluster_exp=np.array(cluster_exp)
    cluster_mean=[np.mean(cluster_exp[i]) for i in range(len(cluster_exp))]
    return variation(cluster_mean)


if __name__ == "__main__":
    formatMatrix()
