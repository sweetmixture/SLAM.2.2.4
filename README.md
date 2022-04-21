## S L A M
Sp Lone pair integrated Atomistic Model - version 2.4.4:   
Semi-empirical force field code for simulating cations including sp-lone pair density.  

For more details, please contact :   
woong.jee.16@ucl.ac.uk / wldndrb1@gmail.com / scott.woodley@ucl.ac.uk

* * *

### 1. Compilation Guide
Makefile is prepared in ```/src``` and for successfuly compilation, it requires items below.  
 - C compiler (the Intel version recommended - especially for the pre-calculated integrals ```*.c``` in ```/src/sp_cluster_integral_lib```).  
 - MPI library (current implementaion is partially parallelised, especially to deal with large number of integral calculations).  
<https://www.open-mpi.org>.  
 - GNU Scientific Library (GSL) for BLAS support.   
<https://www.gnu.org/software/gsl/>.   

For compilation, go to ```/src``` and type ```Make``` will get you the version of executable ```slam.240122.mpi.x``` in the same directory.  
(In the Makefile, please make sure to set appropriate environment variables depend on your system specification)      
* * *

### 2. Standard SLAM Input Files
Before the first SLAM run, contents below explains the basic I/O of the software.
First about the required standard SLAM input files (mandatory) and about their formats.

To familiarise structure of the input and output files, we will take ``` /PbF2_Example1 ```.  
(who cannot really wait, there is a quick guide in ```/PbF2_Example1/ExampleOverview1```)

##### 2.1 List of Mandatory Inputs ☝️
* Basis Set: ```/PbF2_Example1/sp_cluster_parameter_src```.  
SLAM employs a numerical basis set, which describes the molecular orbitals of lone pair electron densities.  
Thus, we have prepared the normalised basis functions. In the directory above, you can find three regular text files,   
```sp_cluster_knot.txt``` and ```sp_cluster_radial_s/p.txt``` where the former and latter contain the knots and the spline functions.  
Those data were taken from the FHI-AIMS code (ab initio software package, <https://fhi-aims.org/>) and re-processed.  
The provided basis set is for Pb(II) cation, and we also provide the sets for other species, Sn(II) and Bi(III), stored in ```/species```
* Geometry : ```/PbF2_Example1/geo.txt```.   
```
#comment
20 5
  F  c    2.341841    2.092431    2.542214
  F  s    2.341841    2.092431    2.542214
  F  c    0.606138    2.515924    4.641656
  F  s    0.606138    2.515924    4.641656
  F  c    0.570914    3.933878    2.286673
  F  s    0.570914    3.933878    2.286673
  F  c    3.500036    3.618994    4.946581
  F  s    3.500036    3.618994    4.946581
  F  c    4.402269    2.779899    1.173653
  F  s    4.402269    2.779899    1.173653
  F  c    1.611292   -0.372446    2.579035
  F  s    1.611292   -0.372446    2.579035
  F  c    2.791786    1.208199   -0.083704
  F  s    2.791786    1.208199   -0.083704
  F  c    3.127108    4.713294    2.528958
  F  s    3.127108    4.713294    2.528958
  F  c    0.273948    1.965956    0.648266
  F  s    0.273948    1.965956    0.648266
  F  c    4.815906    1.248764    3.418853
  F  s    4.815906    1.248764    3.418853
 Pb       2.123537    3.290852    0.550270
 Pb       3.618500    0.274032    1.763194
 Pb       4.942059    3.459499    3.239038
 Pb       0.129201    1.323455    2.777418
 Pb       1.561445    4.477080    4.221179

```
1st line : left as space for comment.   
2nd line :```<var1>   <var2>``` :   
```var1``` specifiy the number of classical species - if the shell model is used, sum of both core(c) and shell(s).    
```var2``` specifiy the number of lone pair speices (cations).  
3rd line ~  
classical species : ```<atom_name> <type> <x> <y> <z> (optional)<_fix_>```  // (optional) fixes the position of species.  
lone-pair species : ```<atom_name> <x> <y> <z> (optional)<_fix_>```  
   
* Potential Parameters : ```/PbF2_Example1/sp_cluster_species.txt```   



