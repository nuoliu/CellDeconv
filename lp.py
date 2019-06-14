"""
Solves linear programming for the deconvolution of RNA-seq expressions

"""
from pulp import *
import pandas as pd
import numpy as np

class Node:
    """
    use a node structure to represent a cell type
    """
    def __init__(self,weights,expr):

        self.left = None
        self.right = None
        self.expr=expr   #keeps the expression vector for this cell type, this an array
        self.weights = weights #keeps an array of weights assigned to this cell type in each sample
                            #self.weights[i] gives the weight assigned by sample with index i
        self.splittable=True
    def getExpression(self):
        return self.expr
    def getWeights(self):
        return self.weights
    def getLeftChild(self):
        return self.left
    def getRightChild(self):
        return self.right

    def isSplittable(self):
        return self.splittable
    def notToSplit(self):
        self.splittable=False

    def __eq__(self, other):
        return ((self.left==other.left) and (self.right==other.right) and np.array_equal(self.weights,other.weights)and np.array_equal(self.expr,other.expr) 
                and (self.splittable==other.splittable))
    def __ne__(self,other):
        return not (self==other)


class CellTree:
    def __init__(self, num_gene, num_samples):
        self.root=Node(None,None)
        self.num_gene=num_gene
        self.num_samples=num_samples
        self.num_leaves=0
        self.leaves=[]  #a list of nodes at the leaf level
        self.addLeaf(self.root)
        self.depth=0
        self.objective=float('inf')    #records the objective achieved by this tree 
    def getNumGenes(self):
        return self.num_gene
    def getObjective(self):
        return self.objective
    def setObjective(self,newObj):
        self.objective=newObj
        

    def splitLeaf(self,leaf,leftChild_weight, leftChild_expr,rightChild_weight,rightChild_expr):
        self.leaves.remove(leaf)
        self.num_leaves=self.num_leaves-1
        leftChild=Node(leftChild_weight,leftChild_expr)
        rightChild=Node(rightChild_weight,rightChild_expr)
        leaf.left=leftChild
        leaf.right=rightChild
        self.addLeaf(leftChild)
        self.addLeaf(rightChild)
        self.depth+=1
        
    
    def addLeaf(self,leaf):
        self.leaves.append(leaf)
        self.num_leaves+=1
    
    def getLeaves(self):
        return self.leaves
    
    def getNumleaves(self):
        return self.num_leaves

    def calculateOffset(self, leaf):
        """
        helper function to calculate the "offset" data matrix if we were to 
        split on the input leaf, i.e the amount of values fixed given the weights of all
        other cell types are fixed and their expression vectors found

        Return
        ---------
        offset:             numpy array of dimension (#gene, #sample)
        """
        offset=np.zeros((self.num_samples,self.num_gene))
        if leaf.getExpression() is not None:
            for i in range(self.num_samples):
                #for each sample, calcuate the weighted sum of all cell type expressions
                for other_leaf in self.leaves:
                    if other_leaf !=leaf:
                        offset[i]=offset[i]+other_leaf.getWeights()[i]*other_leaf.getExpression()
        return offset

            
                





def LP_solver(weight_combo,num_celltype,data,sample_list,num_gene,leaf_to_split, tree):
    if leaf_to_split.getWeights() is None:
        root=True
    else:
        root=False
    offset=tree.calculateOffset(leaf_to_split)


    decon_problem = LpProblem("celltype_Deconvolution", LpMinimize)
    num_sample=len(sample_list)

    expression_vars={}
    for i in range(num_celltype):
        #for each cell type
        #construct a dictionary of variables for the mean value of each gene expression
        #called using expression_vars[celltype_ind][gene_ind]
        expression_vars[i]=LpVariable.dict("cell_type_"+str(i),list(range(num_gene)), lowBound=0, cat="Continuous")

    error_vars={}
    all_errors=[]
    for i in range(num_sample):
        error_vars[i]=LpVariable.dict("error_"+str(i),list(range(num_gene)),lowBound=0,cat="Continuous")
        all_errors.extend(list(error_vars[i].values()))

    #add objective function
    decon_problem+=lpSum(all_errors),"total error"
    #----add constraints----------------
    for i in range(num_sample):
        for j in range(num_gene):
            #calculate the weights for two cell types
            if root:
                weight_to_split=1.0
            else:
                weight_to_split=leaf_to_split.getWeights()[i]
            weight_cell1=weight_combo[i]*weight_to_split
            weight_cell2=(1-weight_combo[i])*weight_to_split

            decon_problem+=expression_vars[0][j]*weight_cell1+expression_vars[1][j]*weight_cell2+offset[i][j]-data[j][i]<=error_vars[i][j]
            decon_problem+=expression_vars[0][j]*weight_cell1+expression_vars[1][j]*weight_cell2+offset[i][j]-data[j][i]>=-error_vars[i][j]
    decon_problem.solve()
    # print("Status of the problem :\n")
    # print(LpStatus[decon_problem.status])
    cell_1_exp=[]
    cell_2_exp=[]
    for i in range(num_gene):
        cell_1_exp.append(expression_vars[0][i].varValue)
        cell_2_exp.append(expression_vars[1][i].varValue)


    return value(decon_problem.objective),np.array(cell_1_exp), np.array(cell_2_exp)




def gridSearchCellTypes(num_celltype=2):
    data=pd.read_excel("data/tumor_filtered_genes.xlsx")
    #drop TAN samples
    data.drop(list(data.filter(regex = 'TAN')), axis = 1, inplace = True)
    #get rid of the annotation and stats columns
    gene_list=data.loc[:,"geneSymbol"].values
    data=data.iloc[:,4:-2]
    tumor_1=data.filter(regex="CG118")
    tumor_2=data.filter(regex="CG163")
    tumor_3=data.filter(regex="CG565")
    # oneTumor("CG118",tumor_1, gene_list)
    # oneTumor("CG163",tumor_2,gene_list)
    oneTumor("CG565",tumor_3,gene_list)



def oneTumor(tumor_name, data,gene_list,num_celltype=2):
    sample_list=data.columns.values
    data=data.values
    num_gene=data.shape[0]
    num_samples=data.shape[1]
    weight_combos=generate_TwoComponentWeights(num_samples)

    #construct a cell tree for the samples
    tree=CellTree(num_gene,num_samples)

    
   
    ####Split the cell types tree#####################################
    while True:
        counter=0
        for leaf in tree.getLeaves():
            if leaf.isSplittable():
                #Try split on this leaf
                best_objective=float("inf")
                best_combo=None
                best_solution=None
                #try all different weight assignments 
                for index,weight_combo in enumerate(weight_combos):
                    # print("Trying the %dth combination of weights"%(index+1))
                    objective, cell_1_exp, cell_2_exp=LP_solver(weight_combo,num_celltype,data,sample_list,num_gene,leaf, tree)
                    if objective<best_objective:
                        best_objective=objective
                        best_combo=weight_combo
                        best_solution=(cell_1_exp,cell_2_exp)
                #whether to use this splitting
                current_obj=tree.getObjective()
                if (current_obj==float('inf')) or ((best_objective<current_obj) and ((current_obj-best_objective)/float(current_obj)>0.05)):
                    #if there is more than 10% decrease in error
                    #adopt this split
                    tree.setObjective(best_objective)
                    leftChild_expr=best_solution[0]
                    rightChild_expr=best_solution[1]
                    total_weights=leaf.getWeights()
                    if total_weights is None: 
                        total_weights=np.array([1.0 for i in range(num_samples)])
                    leftChild_weight=np.multiply(total_weights,np.array(best_combo))
                    complement=np.array([1.0-x for x in best_combo])
                    rightChild_weight=np.multiply(total_weights,complement)
                    
                    tree.splitLeaf(leaf,leftChild_weight,leftChild_expr,rightChild_weight,rightChild_expr)
                else:
                    #if the split is does not improve the fitting
                    leaf.notToSplit()  

            else:
                counter+=1
        if counter==tree.getNumleaves():
            break


    ###finished splitting the cell types#####################


    printCellTypeExpressions(tumor_name, tree,gene_list)
    printSampleWeights(tumor_name, tree, sample_list)
    

def printCellTypeExpressions(tumor_name, tree,gene_list):
    celltypes=tree.getLeaves()
    all_exp=np.zeros((len(celltypes),tree.getNumGenes()))
    for index,celltype in enumerate(celltypes):
        all_exp[index]=celltype.getExpression()
    all_exp=all_exp.transpose()

    vecs=pd.DataFrame(all_exp, index=gene_list,columns=["cell type "+str(i+1) for i in range(len(celltypes))])
    with pd.ExcelWriter('data/results/'+tumor_name+'_inferred_cellTyle_expressions.xlsx') as writer:
        vecs.to_excel(writer, index=False)

def printSampleWeights(tumor_name, tree, sample_list):
    celltypes=tree.getLeaves()
    sample_weights=np.zeros((len(celltypes),len(sample_list)))
    for index,celltype in enumerate(celltypes):
        sample_weights[index]=celltype.getWeights()

    df=pd.DataFrame(sample_weights,index=["cell type "+str(i+1) for i in range(len(celltypes))], columns=sample_list)
    with pd.ExcelWriter('data/results/'+tumor_name+'_inferred_sample_proportions.xlsx') as writer:
        df.to_excel(writer, index=False)






def generate_TwoComponentWeights(num_samples):
    l=[[float(1/3)],[float(2/3)]]
    for i in range(num_samples-1):
        l_new=[]
        for element in l:
           
            new_1=element+[float(1/3)]
            l_new.append(new_1)
     
            new_2=element+[float(2/3)]
            l_new.append(new_2)
   
        l=l_new
    return l




def main():
    gridSearchCellTypes(2)

if __name__=="__main__":
    main()