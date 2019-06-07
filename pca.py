import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data=pd.read_excel('data/genes_higher_interGroup_diff.xlsx')

data=data.drop(["target_id","biotype","geneSymbol","CV","118_v_163","118_v_565","163_v_565"],axis=1)
sample_label=data.columns
data_t=data.T
x=data_t.values
x=StandardScaler().fit_transform(x)
#targets
targets=[x[:5] for x in sample_label]
y=pd.DataFrame([x[:5] for x in sample_label],columns=["target"])

# pca=PCA(n_components=2)
# principalComponents=pca.fit_transform(x)
# principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2'])

# print(pca.explained_variance_ratio_)
# finalDf = pd.concat([principalDf, y], axis = 1)



# fig = plt.figure(figsize = (8,8))
# ax = fig.add_subplot(1,1,1) 
# ax.set_xlabel('Principal Component 1', fontsize = 15)
# ax.set_ylabel('Principal Component 2', fontsize = 15)
# ax.set_title('2 component PCA', fontsize = 20)
# targets = ['CG118', 'CG163', 'CG565']
# colors = ['r', 'g', 'b']
# for target, color in zip(targets,colors):
#     indicesToKeep = finalDf['target'] == target
#     ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
#                , finalDf.loc[indicesToKeep, 'principal component 2']
#                , c = color
#                , s = 50)
# ax.legend(targets)
# ax.grid()
# plt.savefig("PCA_2.jpg")
pca=PCA(n_components=3)
principalComponents=pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents, columns = ['principal component 1', 'principal component 2','principal component 3'])

print(pca.explained_variance_ratio_)
finalDf = pd.concat([principalDf, y], axis = 1)



fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('PC1', fontsize = 15)
ax.set_ylabel('PC2', fontsize = 15)
ax.set_zlabel('PC3', fontsize = 15)
ax.set_title('3 component PCA', fontsize = 20)
targets = ['CG118', 'CG163', 'CG565']
colors = ['r', 'g', 'b']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['target'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'principal component 1']
               , finalDf.loc[indicesToKeep, 'principal component 2']
               , finalDf.loc[indicesToKeep, 'principal component 3']
               , c = color
               , s = 50)
ax.legend(targets)

plt.savefig("PCA_3.jpg")
plt.show()
