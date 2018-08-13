%% Introduction of Vlasov Poisson equation
% We solve Vlasov equation by geralized way
%
% $$\frac{\partial f}{\partial t}+v\nabla_{\bf x} f+E\cdot\nabla_{\bf v}f = 0$$
%
% The Structure of PDE is as follows:
%
%% MATLAB(R) Code
% 
%	for i = 1:10
%       disp x
%	end
% 
%% LaTeX Markup Example
% <latex>
% \begin{tabular}{|r|r|}
% \hline $n$&$n!$\\ 
% \hline 1&1\\ 2&2\\ 3&6\\ 
% \hline
% \end{tabular}
% </latex>
%
%% MATLAB(R) Code
%   formT1.dim = 1;
%	formT1.type = 'FuncMass';
%	formT1.G = @(x)x;
%	formT2.dim = 1;
%	formT2.type = 'FuncGrad';
%	formT2.G = @(x)1;
% 
%	term1 = {formT1,formT2};
% 
%	formT1.dim = 1;
%	formT1.type = 'FuncMass';
%	formT1.G = @(x)x;
% 
%	formT2.dim = 1;
%	formT2.type = 'FuncGrad';
%	formT2.G = @(x)1;
% 
%	term2 = {formT1,formT2};
% 
% 
%	terms = {term1,term2};
% 
%	pde.terms = terms;