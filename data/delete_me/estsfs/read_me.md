Test data to be used for automated testing with the nf-core pipelines


## Introduction
The folder contain three test data files for the module, est-sfs

1). config-JC.txt : This file consists of three lines of text in the following order:
n_outgroup [1, 2, or 3]
model [0, 1 or 2]
nrandom [0 or positive integer]
model
0 = Jukes-Cantor model
1 = Kimura 2-parameter model
2 = Rate-6 model (see Keightley and Jackson 2018 for details)
For further description of this file, refer to the manual of est-sfs. 

2). seedfile.txt
This text file (example provided) contains a single positive integer. It is overwritten by a new random
number seed at the end of the run.

3). TEST-DATA.TXT
An example data file for three outgroups. For further description of this file, refer to the manual of est-sfs. 
