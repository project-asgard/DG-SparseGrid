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
%%% Set initial conditions based on the HashTable
% <InitialCondition.html Initial Conditions>
%% Section 4. Time-independent Matrices
% <TimeIndependentMatrices.html Algorithm 3---Time-independent matrices>
%% Section 5. Time-dependent Matrices
% <TimeDependentMatrices.html Algorithm 4---Time-dependent matrices>
%% Section 6. Global Vlasov Equation
% <GlobalVlasov.html Algorithm 6---global Vlasov equation>
%% Section 7. Time Stepping Method
% <TimeAdvance.html Algorithm 7---time advance method>
%% Section 8. Poisson/Maxwell Solver
%%% Poisson Solver
% <PoissonSolver.html Poisson solver>
%%% Maxwell Solver
% <MaxwellSolver.html Maxwell solver>


