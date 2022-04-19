## S L A M
Sp Lone pair integrated Atomistic Model - version 2.4.4:   
Semi-empirical force field code for simulating cations including sp-lone pair density.  

For more details, please contact :   
woong.jee.16@ucl.ac.uk / wldndrb1@gmail.com / scott.woodley@ucl.ac.uk

* * *

### Compilation Guide
Makefile is prepared in ```/src``` and for successfuly compilation, it requires items below.  
 - C compiler (the Intel version recommended - especially for the pre-calculated integrals ```*.c``` in ```/src/sp_cluster_integral_lib```).  
 - MPI library (current implementaion is partially parallelised, especially to deal with large number of integral calculations).  
<https://www.open-mpi.org>.  
 - GNU Scientific Library (GSL) for BLAS support.   
<https://www.gnu.org/software/gsl/>.   

* * *
### Basic I/O
To familiarise structure of the input and output files, we will take ``` /PbF2_Example1 ```.  
(who cannot really wait, there is a quick guide in ```/PbF2_Example1/ExampleOverview1```)
