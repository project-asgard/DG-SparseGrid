# Matlab Implementation of  Vlasov-Poisson/Maxwell Equations

Matlab version of 
A sparse-grid, Discontinuous Galerkin code for solving the Vlasov-Poisson/Maxwell 
in the general setting (1x1v,1x2v,3x3v).

> 1x1v is for Vlasov-Poisson Solver

> 1x2v if for Vlasov-Maxwell Solver of streaming Weibel (SW) instability

### Run after cloning

```
> run main_Vlasov.m
```

### Run other tests and Modifying the parameters

**Tests can be modified in**
* ./PDE/Vlasov*.m
*%which defines the domain Lmax, Vmax, and Initial condition f(x,v,t=0)*

**Parameters can be modifies in**
* main_Vlasov.m::
 - TEND *%The ending time for test*
 - Lev  *%The mesh's resolution-->with size h=2^(-Lev)*
 - Deg  *%The degree of polynomial--> Deg=2 means linear element*
 - Solver *%Define the solver as VM or VP*
 
* Vlasov*.m::
 - DimX *% The dimensionality of X*
 - DimV *% The dimensionality of V*
 
![Diagram](Vlasov-Maxwell/Notes/VlasovSolver (1).png)
</a href="https://github.com/MollyRaver/drawio-github/blob/master/VlasovSolver.xml"> Edit
