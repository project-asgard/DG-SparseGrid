# Reference matlab implementation of FK6D

A sparse-grid, Discontinuous Galerkin code for solving the Vlasov-Poisson 
in the 2D setting (1x1v).

### Changes in Runaway-electron problem

The file /PDE/Vlasov_RE.m contains the specific functions (listed at the end) and weak-formulation components (listed in the beginning) for the runaway-electron problem. 

### Run after cloning

> run main_Vlasov.m

### Run other tests and Modifying the parameters

**Tests can be modified in**
* ./PDE/Vlasov*.m
*%which defines the domain Lmax, Vmax, and Initial condition f(x,v,t=0)*

**Parameters can be modifies in**
* main_Vlasov.m::
 - TEND *%The ending time for test*
 - Lev  *%The mesh's resolution-->with size h=2^(-Lev)*
 - Deg  *%The degree of polynomial--> Deg=2 means linear element*
 - IsSlowVersion *%IsSlowVersion=1 means using the slow version; other values means using the fast version*
 
* Note::
 - The Hash Table only contains information from Lev
 - The globa matrix is generated through looping over Hash Table
 - Insert the 1D index mapping to HashTable
 - Slow Version of GlobalMatrix_SG is updated in the slowversion foler
