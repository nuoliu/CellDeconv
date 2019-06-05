import pandas as pd
def constructSet():
    L1000=set()
    f = open("data/L1000.txt", "r")
    f.readline()
    for line in f:
        symbol=line.split()[1]
        L1000.add(symbol)
    return L1000


def selectExpression():
#973 of 978 genes from landmark genes are found in the RNA-seq dataset of 
#Wilms tumor, and there are some replicates (1068 in total)
    L1000=constructSet()
    data=pd.read_excel("data/Wilms_gene_tpms_all_samples.xlsx", sheet_name="gene_tpms_all_samples")
    selected=data.loc[data['geneSymbol'].isin(L1000)]
    with pd.ExcelWriter('data/L1000_rnaseq.xlsx') as writer:
        selected.to_excel(writer, index=False)

if __name__=="__main__":
        selectExpression()