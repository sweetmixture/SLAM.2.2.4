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
* Geometry (example) : ```/PbF2_Example1/geo.txt```.   
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
1st line Format : left as space for comment.   
2nd line Format :```<var1>   <var2>``` :   
```var1``` specifiy the number of classical species - if the shell model is used, sum of both core(c) and shell(s).    
```var2``` specifiy the number of lone pair speices (cations).  
3rd line ~  
classical species Format : ```<atom_name> <type> <x> <y> <z> (optional)``` // (optional) `_fix_` the position of species, if specified.  
lone-pair species Format : ```<atom_name> <x> <y> <z> (optional)```  
   
* Potential Parameters (example) : ```/PbF2_Example1/sp_cluster_species.txt```   
```
QM_POT  2
F       122.   0.089
Pb      422.   0.098
ATOM_TYPE 2
SHELL   F       0.59    -1.59   20.77   0.
SP      Pb      3.64    -1.64   14.09
MM_POT  2
SHELL   F       SHELL   F       1127.7  0.2753   -15.83
SP      Pb      SHELL   F       4494.   0.264   -20.4
gnorm_tol 0.000005
```
_Electronic Short-Range Repulsive Potential_, starts with `QM_POT 2`.  

QM_POT 2 Format: `QM_POT <number_interactions>`.  

```F       122.   0.089``` Format : ```<species> <D> <Rho>```  
```Pb      422.   0.098``` Format : ```<species> <D> <Rho>```  

{D,Rho} parameters for the Born-Mayer potential for (LP_density - species).  

_ATOM Species_, starts with `ATOM_TYPE 2`.  

`SHELL   F       0.59    -1.59   20.77   0.`  
`SP      Pb      3.64    -1.64   14.09`  

Format1 (SM ) ```<TYPE=SHELL> <species> <qc> <qs> <k2> <k4>```  
Format2 (RIM) ```<TYPE=RIM> <species> <qc>```                  
Format3 (SP ) ```<TYPE=SP> <speices> <qc> <q_LP> <L>```  

_Standard Classical Interatomic Potential (Buckingham)_, starts with `MM_POT  2`  

`SHELL   F       SHELL   F       1127.7  0.2753   -15.83`  
`SP      Pb      SHELL   F       4494.   0.264   -20.4`  

Format1 ```<TYPE=SP/SHELL/RIM> <species> <TYPE=SHELL/RIM> <speices> <D> <Rho> <-C>```   

Here, note that `SP` refers to the classical part (the core) of the lone pair ion.  
e.g. `SP      Pb` specifies the core part of the lone pair cation `Pb`.   

_Optionals_  
`gnorm_toal <var>` : `var` specifies the termination tolerance - the gradient norm of system.   
i.e. if total force is less than the tolerance, terminates getometry optimisation.   

Others: `h_matrix`, `DIIS`, `get_spline`, `log_limit`, `single` etc. mostly for the developers use.  
* * *
With the above three input files,   
```
/PbF2_Example1/sp_cluster_parameter_src
/PbF2_Example1/geo.txt
/PbF2_Example1/sp_cluster_species.txt
```
now it is able to run the first SLAM run.  

Being in `/PbF2_Example1`, type standard command for running mpi program with the following format, e.g.
(Assuming the excutable is in `/src` as where it was compiled)  

`mpirun -np 2 ./../src/slam.240122.mpi.x <output.xyz> <mx_cycle> (optional)<output.cube> > <standard_out>`     

where  
` <output.xyz>` is a mandatory argument, which will get you the final geometry in *.xyz format.  
`  <mx_cycle>`  is a mandatory argument, one can put any integer to specify the maximum optimisation cycles.  
`<output.cube>` is an optional argument, if you want to get *.cube for visualised lone pair electron densities.  

Say, the command is,   

`mpirun -np 2 ./../src/slam.240122.mpi.x out.xyz 1000 out.cube > out.txt`     

once the run is successfully done, you will get: `out.xyz`, `out.cube`, `out.txt` and `geo.txt.next`,   
where `geo.txt.next` is auto-generated file in the format of `geo.txt`.   
If the optimisation is not fully done, you can simply detach `.next` and feed the dump into a reoptimisation.  
* * *

### 3. Standard Output File

From the above example run, you got the standard output `out.txt` which mainly contains the input file information and the optimisation log.  

For example, taking the very last bit of the output,
```
#############################################################################

 Optimisation Meets Termination Condition, Final Configuration is

 Geometric Derivatives ( eV / Angstrom )

-----------------------------------------------------------------------------
 Species.        x           y           z          |r| 
-----------------------------------------------------------------------------
  F  c   6.410447e-05  -9.005022e-05  -6.022619e-05   1.258794e-04
  F  s   4.286735e-05  -1.153583e-04  -6.963732e-05   1.414019e-04
  F  c   9.984063e-05   1.504981e-05   7.130339e-06   1.012200e-04
  F  s   8.073172e-05   2.151311e-05   1.072822e-05   8.423490e-05
  F  c  -1.366176e-05  -6.132718e-05   4.500437e-05   7.728557e-05
  F  s  -2.215807e-05  -4.190551e-05   4.631848e-05   6.627559e-05
  F  c  -5.054296e-07   6.137396e-05  -9.826536e-07   6.138391e-05
  F  s   4.411737e-07   7.094628e-05   2.031054e-05   7.379761e-05
  F  c   3.125233e-05   1.232911e-05   1.248129e-05   3.583989e-05
  F  s   5.259241e-05   1.937921e-05  -2.310381e-06   5.609682e-05
  F  c  -7.533365e-05  -5.220613e-05  -2.757667e-05   9.571370e-05
  F  s  -8.144294e-05  -5.290119e-05  -2.032788e-05   9.922052e-05
  F  c   3.076065e-05   3.961294e-05  -5.017047e-05   7.093996e-05
  F  s   1.952114e-05   3.159456e-05  -5.467719e-05   6.609755e-05
  F  c  -3.858122e-05  -3.718281e-05  -1.036753e-05   5.457616e-05
  F  s  -3.185434e-05  -1.003380e-06  -1.143972e-05   3.386108e-05
  F  c   2.086171e-06   2.878359e-06  -6.841809e-06   7.710215e-06
  F  s  -3.658755e-06  -2.213751e-07  -3.757834e-06   5.249458e-06
  F  c  -5.012503e-05   2.813586e-05   4.318081e-05   7.189386e-05
  F  s  -2.865585e-05   1.914151e-05   5.086469e-05   6.143917e-05
 Pb      2.575953e-06   5.331827e-06   8.516344e-06   1.037266e-05
 Pb     -6.362898e-05   3.909600e-05   4.305903e-05   8.620455e-05
 Pb     -1.200420e-05   1.065202e-05   1.909900e-06   1.616212e-05
 Pb     -1.076210e-05   1.832829e-05   3.903802e-06   2.160991e-05
 Pb      5.598340e-06   5.679322e-05   2.490782e-05   6.226726e-05
-----------------------------------------------------------------------------

 Cycle              :   969
 Gnorm  (eV/Angs)   :   0.000005
 Energy (eV)        :  -127.530209889502

-----------------------------------------------------------------------------

 Lone Pair Molecular Orbital ( Lowest Eigenvalue / EigenVector )

-----------------------------------------------------------------------------
 Species.   Energy(eV)       s           px           py          pz 
-----------------------------------------------------------------------------
 Pb          9.995407    0.974709    -0.090511    -0.145059     0.143900
 Pb         10.563493    0.973845    -0.158586     0.141812     0.079791
 Pb          9.988816    0.970425    -0.213635    -0.108894    -0.027877
 Pb         10.478729    0.978768     0.130930     0.131765    -0.086655
 Pb         10.555575    0.973756     0.034813    -0.185520    -0.127164
-----------------------------------------------------------------------------

 CONFIGURATION_XYZ_INFO ( final / number of atoms )   :   25
	15
 SCF DONE -127.530210
  F    2.284119    2.240177    2.613017
  F    0.351014    2.617021    4.421014
  F    0.624877    4.483132    2.061168
  F    3.330727    3.483610    4.962918
  F    4.460330    2.596109    1.209927
  F    2.145289   -0.425670    2.493656
  F    2.564866    1.205264   -0.018664
  F    3.227776    4.677504    2.528649
  F    0.073228    1.945430    1.020482
  F    4.753655    1.074218    3.437955
 Pb    1.868984    3.279462    0.590849
 Pb    3.965614    0.323962    1.458518
 Pb    4.858250    3.308968    3.334200
 Pb    0.420859    0.923284    2.959390
 Pb    1.503952    4.490774    4.097544

 CONFIGURATION_XYZ_SC_INFO
 20	5
  F  c    2.284119    2.240177    2.613017
  F  s    2.211383    2.238315    2.547478
  F  c    0.351014    2.617021    4.421014
  F  s    0.419893    2.632530    4.326647
  F  c    0.624877    4.483132    2.061168
  F  s    0.735081    4.428752    2.112381
  F  c    3.330727    3.483610    4.962918
  F  s    3.323589    3.527211    4.830624
  F  c    4.460330    2.596109    1.209927
  F  s    4.434410    2.551059    1.316943
  F  c    2.145289   -0.425670    2.493656
  F  s    2.156757   -0.312751    2.461177
  F  c    2.564866    1.205264   -0.018664
  F  s    2.615843    1.250397    0.095824
  F  c    3.227776    4.677504    2.528649
  F  s    3.227465    4.589219    2.608181
  F  c    0.073228    1.945430    1.020482
  F  s    0.176926    1.947650    1.113830
  F  c    4.753655    1.074218    3.437955
  F  s    4.722330    1.159076    3.331761
 Pb       1.868984    3.279462    0.590849
 Pb       3.965614    0.323962    1.458518
 Pb       4.858250    3.308968    3.334200
 Pb       0.420859    0.923284    2.959390
 Pb       1.503952    4.490774    4.097544

-----------------------------------------------------------------------------
#############################################################################

 Computation Wtime :    198.750706 s

 Date		   :	2022-04-19 06:38:07

 output xyz 	   :	out.xyz
 output cube	   :	out.cube
 next config	   :	geo.txt.next

#############################################################################
```
The file contains some usefule informations,   
e.g. see below: `Geometric Derivatives ( eV / Angstrom )`, `Lone Pair Molecular Orbital ( Lowest Eigenvalue / EigenVector )`,   
which are used for secondary processes:   
 (i) Dipole moment of the system, see `/utils_slam/slam_dipole_analyser`;   
(ii) Normal mode analysis (requires interfacing with Python/Shell written script,   
source: <https://github.com/sweetmixture/slam.vibration_normal_mode.support.git>.  
