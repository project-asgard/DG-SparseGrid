# Matlab Implementation of 2D Vlasov-Poisson Equations

Matlab version of 
A sparse-grid, Discontinuous Galerkin code for solving the Maxwell Equation 
in the 3D setting (3x). Here E=[Ex(x,y,z),Ey(x,y,z),Ez(x,y,z)] and B=[Bx(x,y,z),By(x,y,z),Bz(x,y,z)]

### Run after cloning

> run main_Maxwell.m

### Run other tests and Modifying the parameters

**Tests can be modified in**
* ./PDE/Maxwell*.m
*%which defines the Initial condition E0,B0, and epsilon, mu*

**Parameters can be modifies in**
* main_Maxwell.m::
 - Lev  *%The mesh's resolution-->with size h=2^(-Lev)*
 - Deg  *%The degree of polynomial--> Deg=2 means linear element*
 - pde  *%The test case*
 - MaxT *%time advance steps*
 - dt   *%time stepping*