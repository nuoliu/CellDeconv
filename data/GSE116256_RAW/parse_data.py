import pandas as pd 

all_sample=pd.read_excel("AML328-D113_doubleMalignant.xlsx")
cells=all_sample["Cell"].values
all_exp=pd.read_excel("AML328-D113_expression.xlsx")
cells=cells.tolist()+["Gene"]
all_exp=all_exp.loc[:,cells]
with pd.ExcelWriter("AML328-D113_expression.xlsx") as writer:
    all_exp.to_excel(writer, index=False)
