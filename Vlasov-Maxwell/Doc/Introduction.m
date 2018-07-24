%% Sparse Grids Discontinous Galerkin Methods for Kinetic Equations
%
% SG-DG is a MATLAB software package for solving Vlasov Maxwell/Poisson 
% equations:
%
% $$\frac{\partial f}{\partial t}+v\nabla_{x} f+F\cdot\nabla_{v}f = 0$$
% 
% where $f({\bf x},{\bf v},t)$, and $F={\bf E}+{\bf v}\times{\bf B}$
% (Maxwell Equation) or $F={\bf E}$ (Poisson Equation).
%
% This code is based on the general dimensionality, for example
% 1X1V, 1X2V, and 3X3V systems.
% 


%% Section 1. Introduction
% <NotesVM.html Introduction of Vlasov Maxwell/Poisson equations>
%% Section 2. Numerical Schemes
% <WeakForms.html Weak Forms for VM and VP equations>
%% Section 3. HashTable
% <HashTable.html Algorithm 2---HashTable>
%
% * If we perform the Maxwell or Poisson solver also by SG methods, 
% we will need two sets of Hash. One (Hash1) for $f({\bf x},{\bf v})$ and one (Hash2)
% only contains $x_1,x_2,x_3$ for Maxwell or Poisson equations. 
% However, there will be noncomforming information from these two Hash 
% tables: some grids in Hash2 will not be in Hash1. Not sure whether this
% mismatch will cause trouble or not. 
% * The other way is to generate one Hash table for $f({\bf x},{\bf v})$ but perform
% full grids computing for Poisson or Maxwell solver. This method needs to
% solve a bigger matrix for Poisson or Maxwell. But since we can use
% regular FEM methods and derive ${\bf E}$ and/or ${\bf B}$ only on each element cell,
% it will simplify the calculation of $({\bf E}\nabla_v f,w)$ and/or 
% $(({\bf v}\times {\bf B})\nabla_v f,w)$.
% * We will begin with the full grids calculation for ${\bf E}$ and/or
% ${\bf B}$.
%%% Set initial conditions based on the HashTable
% <InitialConditions.html Initial Conditions>
%
% * We compute the initial conditions for $f({\bf x},{\bf v},0)$ on sparse
% grids. The initial conditions for ${\bf E}$ and/or ${\bf B}$ are for
% scaling functions on the full grids.
%% Section 4. Time-independent Matrices
% <TimeIndependentMatrices.html Algorithm 3---Time-independent matrices>
%
% * The time independet matrices will be generated, including
% ${\bf Grad},\ {\bf NGrad},\ {\bf PGrad},\ {\bf vMassV}$.
%% Section 5. Time-dependent Matrices
% <TimeDependentMatrices.html Algorithm 4---Time-dependent matrices>
% 
% * Compute the time dependent matrix when $E,B$ are known
% ${\bf EMassX,\ BMassX}$
%% Section 6. Global Vlasov Equation
% <GlobalVlasov.html Algorithm 6---global Vlasov equation>
%
% * Assemble the global matrix for Vlasov equations:
% This code works for general dimensionality. For low dimension, some
% parameters $C_{x_1},C_{x_2},C_{x_3},C_{v_1},C_{v_2},C_{v_3}$ are assigned
% to be zeros, and the same assebling techniques
% can be performed.
%% Section 7. Time Stepping Method
% <TimeAdvance.html Algorithm 7---time advance method>
%
% * In this test, we only considered the explicit 3-rd order Runge-Kutta 
% methods. Future work may include the IMEX or implicit time advance
% methods.
%% Section 8. Poisson/Maxwell Solver
%%% Poisson Solver
% <PoissonSolver.html Poisson solver>
%
% * Poisson equations are solved for each time step. When the solution $f_h^{n}$
% is derived, we shall compute ${\bf E}_h^{n+1}=\nabla\phi({\bf x})$.
% * Note:: This calculation is based on the full grids of scaling
% functions.
%%% Maxwell Solver
% <MaxwellSolver.html Maxwell solver>
%
% * Maxwell equations are solved for each time step. When the solution $f_h^{n}$
% is derived, we shall do time involution for ${\bf E}_h^{n+1},{\bf B}_h^{n+1}$
% from the solutions of last time step ${\bf E}_h^{n},{\bf B}_h^{n}$.
% * Note:: This calculation is based on the full grids of scaling
% functions.


