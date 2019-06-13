"""
Solves linear programming for the deconvolution of RNA-seq expressions

"""
from pulp import *
import pandas as pd
import numpy as np

def LP_solver(weight_combo,num_celltype,data,sample_list,num_gene):
    #assume weight_combo is a 1d list with the proportion of cell-type 1
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
            weight_cell1=weight_combo[i]
            weight_cell2=1-weight_cell1
            decon_problem+=expression_vars[0][j]*weight_cell1+expression_vars[1][j]*weight_cell2-data[j][i]<=error_vars[i][j]
            decon_problem+=expression_vars[0][j]*weight_cell1+expression_vars[1][j]*weight_cell2-data[j][i]>=-error_vars[i][j]
    decon_problem.solve()
    print("Status of the problem :\n")
    print(LpStatus[decon_problem.status])
    cell_1_exp=[]
    cell_2_exp=[]
    for i in range(num_gene):
        cell_1_exp.append(expression_vars[0][i].varValue)
        cell_2_exp.append(expression_vars[1][i].varValue)


    return value(decon_problem.objective),cell_1_exp, cell_2_exp

def LP_solver_split(left,other_cell_exp,data,current_splitting,weight_combo,num_celltype,sample_list,num_gene):
    """
    LP solver to split one already established cell type component into two
    """
    #assume weight_combo is a 1d list with the proportion of cell-type 1
    problem = LpProblem("celltype_splitting", LpMinimize)
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


    #add objective function-----------------
    problem+=lpSum(all_errors),"total error"
    #----add constraints----------------
    
    for i in range(num_sample):
        for j in range(num_gene):
            #get the total weight to be split
            #TODO Deal with the situation where the other cell type is already split...
            if isinstance(current_splitting[i],list):
                split_before=True
            else:
                split_before=False
            if left:  #always split left before right
                weight_to_split=current_splitting[i]
            else:   #check if split before to retrieve weight differently
                if split_before:
                    weight_to_split=current_splitting[i][1]
                else:
                    weight_to_split=1-current_splitting[i]
            #get the weight of the cell type not split
            if split_before:
                left_weights=current_splitting[i][0]
                offset=left_weights[0]*other_cell_exp[0][j]+left_weights[1]*other_cell_exp[1][j]
            else:
                fixing_weight=1-weight_to_split
                offset=other_cell_exp[j]*fixing_weight
            #calculate the respective weights of new cell types
            weight_cell1=(weight_combo[i])*weight_to_split
            weight_cell2=(1-weight_combo[i])*weight_to_split
            
            problem+=expression_vars[0][j]*weight_cell1+expression_vars[1][j]*weight_cell2+offset-data[j][i]<=error_vars[i][j]
            problem+=expression_vars[0][j]*weight_cell1+expression_vars[1][j]*weight_cell2+offset-data[j][i]>=-error_vars[i][j]
    problem.solve()
    print("Status of the problem :\n")
    print(LpStatus[problem.status])
    cell_1_exp=[]
    cell_2_exp=[]
    for i in range(num_gene):
        cell_1_exp.append(expression_vars[0][i].varValue)
        cell_2_exp.append(expression_vars[1][i].varValue)


    return value(problem.objective),cell_1_exp, cell_2_exp

def gridSearchCellTypes(num_celltype=2):
    data=pd.read_excel("data/tumor_filtered_genes.xlsx")
    #drop TAN samples
    data.drop(list(data.filter(regex = 'TAN')), axis = 1, inplace = True)
    #get rid of the annotation and stats columns
    data=data.iloc[:,4:-2]
    tumor_1=data.filter(regex="CG118")
    tumor_2=data.filter(regex="CG163")
    tumor_3=data.filter(regex="CG565")
    # oneTumor("CG118",tumor_1)
    oneTumor("CG163",tumor_2)
    oneTumor("CG565",tumor_3)



def oneTumor(tumor_name, data,num_celltype=2):
    sample_list=data.columns.values
    data=data.values
    num_gene=data.shape[0]
    num_samples=data.shape[1]
    weight_combos=generate_TwoComponentWeights(num_samples)
    best_objective=float("inf")
    best_combo=None
    best_solution=None
    for index,weight_combo in enumerate(weight_combos):
        print("Trying the %dth combination of weights"%(index+1))
        objective, cell_1_exp, cell_2_exp=LP_solver(weight_combo,num_celltype,data,sample_list,num_gene)
        if objective<best_objective:
            best_objective=objective
            best_combo=weight_combo
            best_solution=(cell_1_exp,cell_2_exp)
    
    print("The best objective is: %8.3f"% best_objective)
    best_cell_right=best_solution[1]
    best_cell_left=best_solution[0]
    #try splitting on the left cell type------------------
    split_left=split_current(True,best_cell_right,data,best_combo,weight_combos,sample_list, num_gene)
    
    if split_left[0]<best_objective:
        print("YAYYY")
        print("The best objective via splitting on left is: %8.3f"% split_left[0])

        #Update the best objective
        best_objective=split_left[0]

        #update cell type expression vector
        #put the newly found vectors on the left, keep the old right vector, in the form  ((l1,l2),l3)
        best_solution=(split_left[2],best_cell_right)

        #update splitting weights
        best_splitting=split_left[1]   #the best splitting within the firt component
        #a list of lists, each sublist correspond to a sample, has the form [[weight1,weight2],weight3]
        best_combo=[ [ [best_splitting[i]*best_combo[i],(1-best_splitting[i])*best_combo[i]] ,1-best_combo[i]] for i in range(num_samples)]
    
    # try splitting on the right cell type-----------------
    best_cell_left=best_solution[0]    
    split_right=split_current(False, best_cell_left,data, best_combo, weight_combos, sample_list, num_gene)

    if split_right[0]<best_objective:
        print("Yayy")
        print("The best objective via splitting on left is: %8.3f"% split_right[0])
         #Update the best objective
        best_objective=split_right[0]

        #update cell type expression vector
        #put the newly found vectors on the right, keep the old right vector, in the form  ((l1,l2),(l3,l4)) or (l1,(l2,l3))
        best_solution=(best_cell_left,split_right[2])

        #update splitting weights
        best_splitting=split_right[1]   #the best splitting within the second component
        #a list of lists, each sublist correspond to a sample, has the form [[weight1,weight2],weight3] or [[w1,w2],[w3,w4]]
        best_combo=[ [ best_combo[i][0] ,[ best_combo[i][1]*best_splitting[i],  best_combo[i][1]*(1-best_splitting[i])]] for i in range(num_samples)]
    

    printCellTypeExpressions(tumor_name, best_solution)
    printSampleWeights(tumor_name, best_combo, sample_list)
    

def printCellTypeExpressions(tumor_name, best_solution):
    def flatten(vecs):
        copy=vecs
        if isinstance(copy,list):
            return [copy]
        else:
            first=flatten(vecs[0])
            first.extend(flatten(vecs[1]))
            return first
    flattened_vecs=flatten(best_solution)
    vec_array=np.array(flattened_vecs)
    vec_array=vec_array.transpose()
    vecs=pd.DataFrame(vec_array, columns=["cell type "+str(i+1) for i in range(vec_array.shape[1])])
    with pd.ExcelWriter('data/results/'+tumor_name+'_inferred_cellTyle_expressions.xlsx') as writer:
        vecs.to_excel(writer, index=False)

def printSampleWeights(tumor_name, best_combo, sample_list):
    def flatten_sample_weight(weights):
        copy=weights
        if not isinstance(copy,list):
            return [copy]
        else:
            first=flatten_sample_weight(copy[0])
            first.extend(flatten_sample_weight(copy[1]))
            return first
    
    flattened=np.array([flatten_sample_weight(x) for x in  best_combo])
    df=pd.DataFrame(flattened,index=sample_list, columns=["cell type "+str(i+1) for i in range(flattened.shape[1])])
    with pd.ExcelWriter('data/results/'+tumor_name+'_inferred_sample_proportions.xlsx') as writer:
        df.to_excel(writer, index=False)




    
def split_current(left,other_cell_exp,data, current_splitting, weight_combos, sample_list,num_gene):
    """
    other_cell_exp: the expressino vector(list )of the cell type not to be split
    current_splitting: 1D list with the best splitting/weights for 2 cell types
    weight_combos: all the permutations of weights to explore for splitting cell type 1 (left node)
    """

    best_objective=float("inf")
    best_combo=None
    best_solution=None
    for index,weight_combo in enumerate(weight_combos):
        print("Splitting: Trying the %dth combination of weights"%(index+1))
        objective, cell_1_exp, cell_2_exp=LP_solver_split(left,other_cell_exp,data,current_splitting,weight_combo,2,sample_list,num_gene)
        if objective<best_objective:
            best_objective=objective
            best_combo=weight_combo
            best_solution=(cell_1_exp,cell_2_exp)
    return (best_objective, best_combo, best_solution)





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