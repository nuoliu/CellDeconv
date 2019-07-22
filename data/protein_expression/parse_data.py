import pandas as pd
import numpy as np
from scipy.stats import variation
for index in [1]:
    filen="patient"+str(index)+".xlsx"
    data=pd.read_excel(filen)
    num_samples=data.shape[1]-1
    data=data.dropna(thresh=6)
    data=data.fillna(0)

    value=data.values
    value=np.exp(value)
    data.iloc[:,:]=value
    CV=variation(value, axis=1)
    data=data.loc[CV>1.3]
    print(data.shape)
    with pd.ExcelWriter("patient"+str(index)+"filtered.xlsx") as writer:
        data.to_excel(writer, index=True,index_label="geneSymbol")

