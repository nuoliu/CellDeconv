"""
Solves linear programming for the deconvolution of RNA-seq expressions

"""
from pulp import *
import pandas as pd
import numpy as np
import copy 
import sys
import getopt

class Node:
    """
    use a node structure to represent a cell type
    """
    def __init__(self,weights,expr, node_num):

        self.left = None
        self.right = None
        self.node_num=node_num
        self.expr=expr   #keeps the expression vector for this cell type, this an array
        self.weights = weights #keeps an array of weights assigned to this cell type in each sample
                            #self.weights[i] gives the weight assigned by sample with index i
        self.splittable=True
    def getNodeNumber(self):
        return self.node_num
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
        self.root=Node(None,None,1)
        self.num_gene=num_gene
        self.num_samples=num_samples
        self.num_leaves=0
        self.leaves=[]  #a list of nodes at the leaf level
        self.allNodes=[]
        self.addLeaf(self.root)
        self.depth=0
        self.objective=float('inf')    #records the objective achieved by this tree 

    
 
    def getAllNode(self):
        return self.allNodes
    def getRep(self):
        return self.getRep
    def getNumGenes(self):
        return self.num_gene
    def getObjective(self):
        return self.objective
    def setObjective(self,newObj):
        self.objective=newObj
        

    def splitLeaf(self,leaf,leftChild_weight, leftChild_expr,rightChild_weight,rightChild_expr):
        self.leaves.remove(leaf)
        self.num_leaves=self.num_leaves-1
        num_Nodes=len(self.getAllNode())
        leftChild=Node(leftChild_weight,leftChild_expr,num_Nodes+1)
        rightChild=Node(rightChild_weight,rightChild_expr, num_Nodes+2)
        leaf.left=leftChild
        
        leaf.right=rightChild
        # if leaf.getLeftChild() is not None:
        #     print("Node %d has children Node %d and Node %d"%(leaf.getNodeNumber(),num_Nodes+1,num_Nodes+2))
        self.addLeaf(leftChild)
        self.addLeaf(rightChild)
        self.depth+=1
        
    
    def addLeaf(self,leaf):
        self.leaves.append(leaf)
        self.allNodes.append(leaf)
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

            
                





def LP_solver(weight_combo,data,sample_list,num_gene,leaf_to_split, tree,num_celltype=2):

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




def gridSearchCellTypes(input_fn, out_dir, fine_grid=True, alpha=1.0):
    # data=pd.read_excel("data/tumor_filtered_genes.xlsx")
    # #drop TAN samples
    # data.drop(list(data.filter(regex = 'TAN')), axis = 1, inplace = True)
    # #get rid of the annotation and stats columns
    # gene_list=data.loc[:,"geneSymbol"].values
    # data=data.iloc[:,4:-2]
    # tumor_1=data.filter(regex="CG118")
    # tumor_2=data.filter(regex="CG163")
    # tumor_3=data.filter(regex="CG565")
    # # oneTumor("CG118",tumor_1,gene_list, True)
    # # oneTumor("CG163",tumor_2,gene_list,True)
    # oneTumor("CG565",tumor_3,gene_list, True)

    data=pd.read_excel(input_fn)
    gene_list=data.loc[:,"geneSymbol"].values
    data=data.iloc[:,1:]
    project_name=input_fn.split('/')[-1].split('.')[0]
    oneTumor(project_name,data,gene_list, out_dir, fine_grid, alpha)


def oneTumor(tumor_name, data,gene_list,out_dir, fine_grid, alpha):
    """
    The actual function that does the grid search and cell type splits
    """
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
        leaves=copy.copy(tree.getLeaves())
        #allows only four leaves
        
        for leaf in leaves:
            if leaf.isSplittable():
                #Try split on this leaf
                best_objective=float("inf")
                best_combo=None
                best_solution=None
                #try all different weight assignments 
                for index,weight_combo in enumerate(weight_combos):
                    # print("Trying the %dth combination of weights"%(index+1))
                    objective, cell_1_exp, cell_2_exp=LP_solver(weight_combo,data,sample_list,num_gene,leaf, tree)
                    if objective<best_objective:
                        best_objective=objective
                        best_combo=weight_combo
                        best_solution=(cell_1_exp,cell_2_exp)
                #whether to use this splitting
                current_obj=tree.getObjective()
                cells_made=tree.getNumleaves()
                best_objective=scaledError(best_objective,cells_made, num_samples, alpha)
                if  best_objective<current_obj:   #compared with penalized error when deciding on splitting
                    #if there is more than 10% decrease in error
                    #adopt this split
                    if fine_grid:
                        #TODO fine tune the splitting we found that is good
                        for i in range(num_samples):
                 
                            #try fine tune the weight in each sample
                            old_combo=best_combo
                            new_combo=refinedGrid(best_combo,i)
                            for combo in new_combo:
                                #try making the weight smaller or larger
                                # print("Trying rid for %dth sample"%(i+1))
                                objective, cell_1_exp, cell_2_exp=LP_solver(combo,data,sample_list,num_gene,leaf, tree)
                                scaled_new_error=scaledError(objective, cells_made,num_samples,alpha)
                                if scaled_new_error<0.9*best_objective:
                                    best_objective=scaled_new_error
                                    best_combo=combo
                                    best_solution=(cell_1_exp,cell_2_exp)
                            if best_combo!=old_combo:
                                #refine one more time
                                new_combo=refinedGrid(best_combo,i)
                                for combo in new_combo:
                                #try making the weight smaller or larger
                                    # print("Trying rid for %dth sample"%(i+1))
                                    objective, cell_1_exp, cell_2_exp=LP_solver(combo,data,sample_list,num_gene,leaf, tree)
                                    scaled_new_error=scaledError(objective, cells_made,num_samples,alpha)
                                    if scaled_new_error<0.9*best_objective:
                                        best_objective=objective
                                        best_combo=combo
                                        best_solution=(cell_1_exp,cell_2_exp)
                  

                    # print("The old objective is %8.3f"%tree.getObjective())
                    # print("The new objective is %8.3f"%best_objective)
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
                    if best_objective==0.0:
                        break
                else:
                    #if the split is does not improve the fitting
                    leaf.notToSplit()  

            else:
                counter+=1
        if counter==tree.getNumleaves() or best_objective==0.0:
            break


    ###finished splitting the cell types#####################

    print("The resulting objective is %8.3f"%tree.getObjective())
    printTree(tree, gene_list, tumor_name, out_dir)
    printCellTypeExpressions(tumor_name, tree,gene_list, out_dir)
    printSampleWeights(tumor_name, tree, sample_list,out_dir)
    
# def isImprove(newObj, oldObj, cells_made, num_samples):
#     """
#     see if the new objective count as an improvement from the old one
#     """
    # return (newObj<oldObj) and ((oldObj-newObj)/float(oldObj)>thr)
def scaledError(newObj, cells_made, num_samples,alpha):
    """
    scales up the total error from linear programming by a penalty factor
    error=error/(1-cells_made/(alpha*num_samples)) 
    """
    denominator=max(0.0001,float(1.0-cells_made/float(alpha*num_samples)))
    return newObj/denominator

def printCellTypeExpressions(tumor_name, tree,gene_list,out_dir):
    celltypes=tree.getLeaves()
    all_exp=np.zeros((len(celltypes),tree.getNumGenes()))
    for index,celltype in enumerate(celltypes):
        all_exp[index]=celltype.getExpression()
    all_exp=all_exp.transpose()

    vecs=pd.DataFrame(all_exp,columns=["cell type "+str(i+1) for i in range(len(celltypes))])
    vecs.insert(0,"geneSymbol", gene_list)
   
    filename=out_dir+tumor_name+'_LP_cellTyle_expressions.xlsx'
    
    with pd.ExcelWriter(filename) as writer:
        vecs.to_excel(writer, index=False)

def printSampleWeights(tumor_name, tree, sample_list,out_dir):
    celltypes=tree.getLeaves()
    sample_weights=np.zeros((len(celltypes),len(sample_list)))
    for index,celltype in enumerate(celltypes):
        sample_weights[index]=celltype.getWeights()

    df=pd.DataFrame(sample_weights, columns=sample_list)
    df.insert(0,'cell type',["cell type "+str(i+1) for i in range(len(celltypes))])
    
    filename=out_dir+tumor_name+'_inferred_proportions.xlsx'
    with pd.ExcelWriter(filename) as writer:
        df.to_excel(writer, index=False)


def generate_TwoComponentWeights(num_samples):
    l=[[float(1/4)],[float(3/4)]]
    for i in range(num_samples-1):
        l_new=[]
        for element in l:
           
            new_1=element+[float(1/4)]
            l_new.append(new_1)
     
            new_2=element+[float(3/4)]
            l_new.append(new_2)
   
        l=l_new
    return l

def refinedGrid(best_combo,sample_ind):
    """
    Input
    ---------------
    best_combo:         a list where each element tells the relative weight 
                        assigned to cell 1 among 2 cell types in each sample
                        example: [0.25, 0.75, 0.75, 0.25, 0.25]

    Output
    ---------------
    refined_combos:     2 modified lists with the element at sample_ind decreased/increased
    """
    def finer_weight(original_weight):
        """
        helper function. Input is a weight, output is a list of two weights around the input
        """
        if original_weight<0.5:
            if original_weight<=0.25:
                return [0.5*original_weight,1.5*original_weight]
            else:
                rest=0.5-original_weight
                return [original_weight-0.5*rest,original_weight+0.5*rest]
        else:
            if original_weight>=0.75:
                rest=1-original_weight
                return [original_weight-0.5*rest,original_weight+0.5*rest]
            else:
                rest=original_weight-0.5
                return [original_weight-0.5*rest,original_weight+0.5*rest]

    #for each sample
    l1=best_combo
    l2=best_combo
    l1[sample_ind]=finer_weight(best_combo[sample_ind])[0]
    l2[sample_ind]=finer_weight(best_combo[sample_ind])[1]
    return (l1,l2)

def printTree(tree, gene_list, tumor_name, out_dir):
    all_nodes=tree.getAllNode()
    num_nodes=len(all_nodes)
    all_exp=[]
    node_names=[]

    for index,node in enumerate(all_nodes):
        it=node.getNodeNumber()
        if node.getExpression() is not None:
            all_exp.append(node.getExpression())
            node_names.append(node.getNodeNumber())
        
    all_exp=np.array(all_exp)

    all_exp=all_exp.transpose()


    vecs=pd.DataFrame(all_exp,columns=node_names)
    vecs.insert(0,"geneSymbol", gene_list)
    
    filename=out_dir+tumor_name+'_expression_allNodes.xlsx'
    with pd.ExcelWriter(filename) as writer:
        vecs.to_excel(writer, index=False)

 
    
   

        

def main(argv):
    inputfn=""
    output_dir=""
    fine_grid=True
    alpha=None
    try:
        opts, args = getopt.getopt(argv,"hi:o:f:a:")
    except getopt.GetoptError:
        print('lp.py -i <input excel file> -o <out_directory> -f <True or False> -a alpha(float)')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('lp.py -i <input excel file> -o <out_directory> -f <True or False> -a alpha(float)')
            sys.exit()
        elif opt=="-i":
            inputfn=arg
        elif opt=="-o":
            output_dir=arg
        elif opt=="-f":
            fine_grid=arg=="True"
        elif opt=="-a":
            alpha=float(arg)
    
  
    gridSearchCellTypes(inputfn,output_dir,fine_grid,alpha)

if __name__=="__main__":
    main(sys.argv[1:])