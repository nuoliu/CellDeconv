import pandas as pd
def get_sample_names():
    f=open("data/tmp2/CCLE/sample_annotation.txt","r")
    f.readline()
    samples=set()
    for line in f:
        info=line.split()
        if "Kidney" in info[2]:
            samples.add(info[0])
    return samples

def select_sample_data():
    samples=get_sample_names()
    data=pd.read_excel("data/tmp2/CCLE/CCLE_data.xlsx")
    col_list=[0]
    for index, name in enumerate(data.columns):
        if name in samples:
            col_list.append(index)
    print(col_list)
    selected_cols=data.iloc[:,col_list]
    selected_cols.rename(columns={"GeneSymbol_DataType":"geneSymbol"},inplace=True)
#     with pd.ExcelWriter('data/kidneyCCL.xlsx') as writer:
#         selected_cols.to_excel(writer, index=False)
    #selecting only data for genes with high variation in wilms tumor
#     high_CV=pd.read_excel("data/genes_highCV.xlsx")
#     genes=high_CV[["geneSymbol"]]
#     print(genes.shape[0])

#     filtered=genes.merge(selected_cols,how="inner",on="geneSymbol")
#     with pd.ExcelWriter('data/kidneyCCL_highCVgenes.xlsx') as writer:
#         filtered.to_excel(writer, index=False)
#     genes=filtered[["geneSymbol"]]
#     high_CV=genes.merge(high_CV,how="inner",on="geneSymbol")
#     #rewrite
#     with pd.ExcelWriter('data/genes_highCV.xlsx') as writer:
#         high_CV.to_excel(writer, index=False)
    data=pd.read_excel("data/L1000_tumor_rnaseq.xlsx")
    data.dropna()
    #deal with duplicates:
    data.drop_duplicates(subset="geneSymbol",inplace=True)
    genes=data[["geneSymbol"]]
    filtered=genes.merge(selected_cols,how="inner",on="geneSymbol")
    with pd.ExcelWriter('data/kidneyCCL_common.xlsx') as writer:         
        filtered.to_excel(writer, index=False)
    genes=filtered[["geneSymbol"]]
    common_mixture=genes.merge(data,how="inner",on="geneSymbol")
    #rewrite
    with pd.ExcelWriter('data/tumor_mixture_common.xlsx') as writer:
        common_mixture.to_excel(writer, index=False)
        

#     genes=set(high_CV.geneSymbol.values)
#     filtered=selected_cols.loc[selected_cols["GeneSymbol_DataType"] in genes]
   




def main():
    select_sample_data()

if __name__=="__main__":
    main()