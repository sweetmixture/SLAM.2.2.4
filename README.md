## SLAM
Sp Lone pair integrated Atomistic Model ( S L A M ) version 2.4.4:   
Semi-empirical force field code for simulating cations including sp-lone pair density.  

For more details, please contact :   
woong.jee.16@ucl.ac.uk / wldndrb1@gmail.com / scott.woodley@ucl.ac.uk

* * *

### Compilation Guide
Makefile is prepared in ```/src``` and for successfuly compilation, it requires items below.  
 - C compiler (the Intel version recommended - especially for the pre-calculated integrals ```*.c``` in ```/src/sp_cluster_integral_lib```).  
 - MPI library (current implementaion is partially parallelised, especially to deal with large number of integral calculations).  
 - GNU Scientific Library for BLAS support.  

* * *
### Basic I/O
To familiarise structure of the input and output files, we will take ``` /PbF2_Example1 ```
