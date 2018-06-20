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

<a href="http://jgraph.github.io/drawio-github/edit-diagram.html?repo=drawio-github&path=diagram.png" target="_blank">Edit</a> | <a href="https://www.draw.io/#Uhttps%3A%2F%2Fjgraph.github.io%2Fdrawio-github%2Fdiagram.png" target="_blank">Edit As New</a>

<a href="http://jgraph.github.io/drawio-github/edit-diagram.html" target="_blank">edit-diagram.html</a> does the I/O with GitHub and uses draw.io in embed mode for diagram editing. The page supports the following URL parameters: user, pass, repo, path, ref and action=open (the Edit link above is an example). Using action=open, links for immediate diagram editing in GitHub can be created (requires user and pass parameters). You can also use files on GitHub as templates in draw.io via the url parameter (see Edit As New above).