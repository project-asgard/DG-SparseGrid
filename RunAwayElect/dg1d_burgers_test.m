function dg1d_burgers_test ( )

%*****************************************************************************80
%
%% DG1D_BURGERS_TEST is a driver for solving the 1D Burgers equation.
%
%  Licensing:
%
%    Permission to use this software for noncommercial
%    research and educational purposes is hereby granted
%    without fee.  Redistribution, sale, or incorporation
%    of this software into a commercial product is prohibited.
%
%    THE AUTHORS OR PUBLISHER DISCLAIMS ANY AND ALL WARRANTIES
%    WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
%    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
%    PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHORS OR
%    THE PUBLISHER BE LIABLE FOR ANY SPECIAL, INDIRECT OR
%    CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
%    RESULTING FROM LOSS OF USE, DATA OR PROFITS.
%
%  Modified:
%
%    17 September 2018
%
%  Author:
%
%    Original version by Jan Hesthaven, Tim Warburton.
%    Some modifications by John Burkardt.
%
%  Reference:
%
%    Jan Hesthaven, Tim Warburton,
%    Nodal Discontinuous Galerkin Methods: 
%    Algorithms, Analysis, and Applications,
%    Springer, 2007,
%    ISBN: 978-0387720654.
%
  addpath ( '../dg1d_burgers' );

%   timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DG1D_BURGERS_TEST:\n' );
  fprintf ( 1, '  MATLAB version.\n' );
  fprintf ( 1, '  Discontinuous Galerkin method for the 1D Burgers equation.\n' );
%
%  Define global variables.
%
  Globals1D;
%
%  Set the polynomial order.
%
  N = 4;
%
%  Generate simple mesh over [-1,+1].
%
  xL = -1.0; 
  xR = 1.0;
  [ Nv, VX, K, EToV ] = MeshGen1D ( xL, xR, 8 );
%
%  Initialize solver and construct grid and metric
%
  StartUp1D;
%
%  Set u at the initial time.
%
  epsilon = 0.1;
  u = - tanh ( ( x + 0.5 ) / ( 2 * epsilon ) ) + 1.0;
%
%  Rearrange the data as vectors.
%
  xv = reshape ( x, (N+1)*8, 1 );
  uv = reshape ( u, (N+1)*8, 1 );
%
%  Plot the solution.
%
  plot ( xv, uv, 'linewidth', 3 );
  grid ( 'on' );
  xlabel ( '<---X--->' )
  ylabel ( '<--U(X,T)-->' )
  title ( 'Burgers Equation Solution at Initial Time' )
%
%  Save the plot in a file.
%
  filename = 'dg1d_burgers_test_initial.png';
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Graphics saved in file "%s"\n', filename );
%
%  Solve the problem.
%
  FinalTime = 1.5;
  [ u ] = Burgers1D ( u, epsilon, xL, xR, FinalTime );
%
%  Rearrange the data as vectors.
%
  xv = reshape ( x, (N+1)*8, 1 );
  uv = reshape ( u, (N+1)*8, 1 );
%
%  Plot the solution.
%
  plot ( xv, uv, 'linewidth', 3 );
  grid ( 'on' );
  xlabel ( '<---X--->' )
  ylabel ( '<--U(X,T)-->' )
  title ( 'Burgers Equation Solution at Final Time' )
%
%  Save the plot in a file.
%
  filename = 'dg1d_burgers_test_final.png';
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Graphics saved in file "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DG1D_BURGERS_TEST:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( )

  rmpath ( '../dg1d_burgers' );

  return
end

