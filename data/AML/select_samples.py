import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
import random
all_sample=pd.read_excel("AML328-D113_doubleMalignant.xlsx")
all_exp=pd.read_excel("AML328-D113_expression.xlsx")
gene_list=all_exp["Gene"].values
cell_types=all_sample["CellType"].unique()
keep_sample=[]
centers=[]
for celltype in cell_types:
    data=all_sample.loc[all_sample["CellType"]==celltype]
    cells=data["Cell"].values.tolist()
    if len(cells)<5:
        keep_sample.extend(cells)
        #compute the average!
        celltype_exp=all_exp.loc[:,cells]
        expression=celltype_exp.values.T
        avg=np.mean(expression,axis=0)
        centers.append(avg)
    else:
        #do k means
        celltype_exp=all_exp.loc[:,cells]
        expression=celltype_exp.values.T
        cell_labels=celltype_exp.columns.values
        kmeans = KMeans(n_clusters=3, random_state=0).fit(expression)
        #select cluster
        unique_clusters, counts_elements = np.unique(kmeans.labels_, return_counts=True)
        #keep only cluster with more than 1/3 of the members
        unique_clusters=[unique_clusters[x] for x in range(3) if counts_elements[x]>len(cells)/3]
        chosen_cluster=random.choice(unique_clusters)
        #pick out the cells that were in this cluster
        indices=np.argwhere(kmeans.labels_==chosen_cluster)
        center=kmeans.cluster_centers_[chosen_cluster]
        for i in indices:  #i is a list with one index
            keep_sample.append(cell_labels[i[0]])
        centers.append(center)
centers=np.array(centers)
centers=centers.T
centers_df=pd.DataFrame(centers, columns=cell_types)
centers_df.insert(0,"geneSymbol",gene_list)
with pd.ExcelWriter("AML328-D113_trueCenters.xlsx") as writer:
    centers_df.to_excel(writer, index=False)


selected_samples=all_sample.loc[all_sample["Cell"].isin(keep_sample)]
with pd.ExcelWriter("AML328-D113_selectedSamples.xlsx") as writer:
    selected_samples.to_excel(writer, index=False)


        




