import pandas as pd
from scipy.stats import variation
from scipy.stats import pearsonr
from collections import *
import numpy as np
# def calculate_intergroup():
#     data=pd.read_excel("data/L1000_tumor_rnaseq.xlsx")
#     nrows=data.shape[0]
#     tumor_118=data[[col for col in data if "118" in col]]
#     tumor_163=data[[col for col in data if "163" in col]]
#     tumor_565=data[[col for col in data if "565" in col]]
#     t_1=[]
#     t_2=[]
#     t_3=[]
#     for i in range(nrows):
#         t,p=stats.ttest_ind(tumor_118.iloc[i],tumor_163.iloc[i],equal_var=False)
#         t_1.append(p)
#         t,p=stats.ttest_ind(tumor_118.iloc[i],tumor_565.iloc[i],equal_var=False)
#         t_2.append(p)
#         t,p=stats.ttest_ind(tumor_163.iloc[i],tumor_565.iloc[i],equal_var=False)
#         t_3.append(p)

#     data["118_v_163"]=t_1
#     data["118_v_565"]=t_2
#     data["163_v_565"]=t_3

#     with pd.ExcelWriter('data/L1000_tumor_with_variation.xlsx') as writer:
#             data.to_excel(writer, index=False)


def select_significant():
    data=pd.read_excel('data/L1000_tumor_rnaseq.xlsx')
    data.dropna()
    #deal with duplicates:
    data.drop_duplicates(subset="geneSymbol",inplace=True)
    #select genes with average expression level high enough
    mean_expression=data["mean_expr"].mean()
    std_expression=data['mean_expr'].std()
    thr=mean_expression-0.2*std_expression
    above_thr_exp=data.loc[data["mean_expr"]>=thr]

    #select genes with >0.3 CV
    higher_variation=above_thr_exp.loc[above_thr_exp["CV"]>0.5]


    #filter out genes which has spearman correlation coeff=1 with other genes
    row_to_keep=list(range(higher_variation.shape[0]))
 
    for i in range(higher_variation.shape[0]):
        
        #first 3 are annotation and last two are stats
        current_gene=higher_variation.iloc[i,3:-2].values
        for j in row_to_keep:
            if j!=i:
                other_gene=higher_variation.iloc[i,3:-2].values
                r,p=stats.spearmanr(current_gene,other_gene)
                if r==1:
                    row_to_keep.remove(i)
                    break

    non_correlated=higher_variation.iloc[row_to_keep]

    #leaves 199 genes
    with pd.ExcelWriter('data/tumor_filtered_genes.xlsx') as writer:
        non_correlated.to_excel(writer, index=False)
                        



    # much_intra_diff=data.loc[(data["118_v_163"]<0.05) &(data["118_v_565"]<0.05) &(data["163_v_565"]<0.05)]
    # high_CV=data.sort_values(by=["CV"],ascending=False)[:100]
    # with pd.ExcelWriter('data/genes_highCV.xlsx') as writer:
    #         high_CV.to_excel(writer, index=False)
    # print(high_CV.iloc[99]["CV"])
    # print(much_intra_diff.shape[0])
#     with pd.ExcelWriter('data/genes_higher_interGroup_diff.xlsx') as writer:
#             much_intra_diff.to_excel(writer, index=False)

def filterWilmsGenes():
    data=pd.read_excel('L1000_tumor_rnaseq.xlsx')
    data.dropna()
    #deal with duplicates:
    data.drop_duplicates(subset="geneSymbol",inplace=True)
    tumor_118=data[["geneSymbol"]+[col for col in data if "118" in col]]
    tumor_163=data[["geneSymbol"]+[col for col in data if "163" in col]]
    tumor_565=data[["geneSymbol"]+[col for col in data if "565" in col]]
    for tumor,name in [(tumor_118,"tumor_118"),(tumor_163,"tumor_163"),(tumor_565,"tumor_565")]:
        expressions=tumor.iloc[:,1:].values
        #filter out genes with max lower than 1
        MAX=np.amax(expressions, axis=1)
        tumor=tumor.loc[MAX>1.0]
        expressions=tumor.iloc[:,1:].values
        CV=variation(expressions, axis=1)
        tumor=tumor.loc[CV>1.0]
        print(name)
        print(tumor.shape)
        with pd.ExcelWriter(name+".xlsx") as writer:
            tumor.to_excel(writer, index=False)


def main():
    filterWilmsGenes()

if __name__=="__main__":
    main()



