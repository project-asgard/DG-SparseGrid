# Matlab Implementation of 2D Vlasov-Poisson Equations

Matlab version of 
A sparse-grid, Discontinuous Galerkin code for solving the Vlasov-Poisson 
in the 2D setting (1x1v).

Run
after cloning

run main_Vlasov.m

Modifying and Run other tests

Parameters can be modified in
./PDE/Vlasov*.m
% which defines the domain Lmax, Vmax, and Initial condition f(x,v,t=0)

main_Vlasov.m::
 TEND % The ending test for test
 Lev  % The mesh's resolution-->with size h=2^(-Lev)
 Deg  % The degree of polynomial--> Deg=2 means linear element