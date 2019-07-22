# CellDeconv
Using a hierachical linear programming approach to breakdown mixture tumor sample expressin profile into multiple cell type expressions.

Author: Nuo (Ivy) Liu, Harvey Mudd College Class of 2020, during the 2019 SMART program at BCM in the lab of Dr. Pavel Sumazin.


## Main programs:
_lp.py_: program that runs the linear programming on input mixture samples, generate three output files: inferred cell type expression (xlsx), all cell type expressions (include internal nodes) (xlsx), and composition/proportion of each inferred cell type in each sample (xlsx). Alpha is a float value around 1 that constrains how many cell types can be produced compared to the number of samples. The higher the alpha, the more number of cell types possible, try testing between 0.4 to 1.6. 

_goodness.py_: If given ground truth (true cell types used to make mixture, etc.) can use this this program to measure the goodness of the inference. Gives a averaged R score between all infered cell types and the true cell types, also geneates baseline R value (using average of mixture sample inputs), and random synthetic cell line R value for comparison.

Usage:
```
python lp_2.py -i PATH_TO_INPUT_MIXTURE_SAMPLES.xlsx -o OUT_DIR/ -f True -a ALPHA
python goodness.py TRUE_CELL_LINE_EXPRESSIONS.xlsx INFERED_CELL_LINE_EXPRESSIONS.xlsx MIXTURE_SAMPLE_EXPRESSIONS.xlsx WORK_DIR/
```

## Directories:
### data
The data used for each experiment is organized in their own folder

### code:
All the python scripts are organized there, some python scripts for preprocessing data are in their own individual data folder if not found here. You can find scripts for generating synthetic mixture data.

### Shell scripts:
All the linux shell scripts use to run tests are organized here, see content to view the requirements and dependencies, the organization of the current directories might not work 

