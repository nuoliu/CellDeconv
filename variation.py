import pandas as pd
from scipy import stats
def calculate_intergroup():
    data=pd.read_excel("data/L1000_tumor_rnaseq.xlsx")
    nrows=data.shape[0]
    tumor_118=data[[col for col in data if "118" in col]]
    tumor_163=data[[col for col in data if "163" in col]]
    tumor_565=data[[col for col in data if "565" in col]]
    t_1=[]
    t_2=[]
    t_3=[]
    for i in range(nrows):
        t,p=stats.ttest_ind(tumor_118.iloc[i],tumor_163.iloc[i],equal_var=False)
        t_1.append(p)
        t,p=stats.ttest_ind(tumor_118.iloc[i],tumor_565.iloc[i],equal_var=False)
        t_2.append(p)
        t,p=stats.ttest_ind(tumor_163.iloc[i],tumor_565.iloc[i],equal_var=False)
        t_3.append(p)

    data["118_v_163"]=t_1
    data["118_v_565"]=t_2
    data["163_v_565"]=t_3

    with pd.ExcelWriter('data/L1000_tumor_with_variation.xlsx') as writer:
            data.to_excel(writer, index=False)


def select_significant():
    data=pd.read_excel('data/L1000_tumor_with_variation.xlsx')
    data.dropna()
    much_intra_diff=data.loc[(data["118_v_163"]<0.05) &(data["118_v_565"]<0.05) &(data["163_v_565"]<0.05)]
    high_CV=data.sort_values(by=["CV"],ascending=False)[:100]
    with pd.ExcelWriter('data/genes_highCV.xlsx') as writer:
            high_CV.to_excel(writer, index=False)
    print(high_CV.iloc[99]["CV"])
    # print(much_intra_diff.shape[0])
    # with pd.ExcelWriter('data/genes_higher_interGroup_diff.xlsx') as writer:
    #         much_intra_diff.to_excel(writer, index=False)
def main():
    select_significant()

if __name__=="__main__":
    main()



