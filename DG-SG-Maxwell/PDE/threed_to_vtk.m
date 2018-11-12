function threed_to_vtk ( x, e_conn, data, output_filename, title )

%% THREED_TO_VTK writes out a 3D FEM scalar solution to a legacy VTK file.
%
%  Discussion:
%
%    The VTK file can be read and displayed by the Paraview program.
%
%  Author:
%
%    Lin Mu
%
%  Parameters:
%
%    Input, real X(NODE_NUM,3), the node coordinates.
%
%    Input, integer E_CONN(ELEMENT_NUM,ELEMENT_ORDER), the nodes that
%    form each element.
%
%    Input, real data(NODE_NUM,1), data at each node.
%
%    Input, string OUTPUT_FILENAME, the name of the output file.
%    By convention, this file should have the extension ".vtu".
%
%    Input, string TITLE, a title for the data.
%
  [ node_num, dim_num ] = size ( x );
  [ element_num,  element_order ] = size ( e_conn );
%
%  Open the output file and write the data to it.
%
  if ( isempty ( output_filename ) )
    output_filename = '3d_fem.vtk';
  end

  output_unit = fopen ( output_filename, 'w' );

  vtk_puv_write ( output_unit, title, node_num, element_num, ...
    element_order, x', e_conn', data);

  fclose ( output_unit ) ;

  fprintf ( 1, '\n' );
  fprintf ( 1, 'THREED_TO_VTK\n' );
  fprintf ( 1, '  Solution data written to "%s"\n', output_filename );

  return
end

function vtk_puv_write ( output_unit, title, node_num, element_num, ...
  element_order, xy, element_node, data )

%% VTK_PUV_WRITE writes scalar data to a VTK file.
%
%  Discussion:
%
%    The data is assumed to have been computed by a finite element program
%    for a 3D geometry which has been meshed using tetrahedral elements 
%
%    Lin Mu
%
%  Parameters:
%
%    Input, integer OUTPUT_UNIT, the output unit.
%
%    Input, string TITLE, a title for the data.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer ELEMENT_ORDER, the order of the elements.
%
%    Input, real XY(3,NODE_NUM), the node coordinates.
%
%    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM), the
%    nodes that make up each element.
%
%    Input, real data(NODE_NUM), the data at each node.
%
%

  fprintf ( output_unit, '# vtk DataFile Version 2.0\n' );
  fprintf ( output_unit, '%s\n', title );
  fprintf ( output_unit, 'ASCII\n' );
  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'DATASET UNSTRUCTURED_GRID\n' );
  fprintf ( output_unit, 'POINTS %d double\n', node_num );

  for node = 1 : node_num
    fprintf ( output_unit, '  %f  %f  %f\n', xy(:,node) );
  end
%
%  Note that CELL_SIZE uses ELEMENT_ORDER+1.
%
%  Note that the 1-based node indices in ELEMENT_NODE must be
%  converted to 0-based indices before being written out.
%
  cell_size = element_num * ( element_order + 1 );

  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'CELLS  %d  %d\n', element_num, cell_size );
  for element = 1 : element_num
    fprintf ( output_unit, '  %d', element_order );
    for order = 1 : element_order
      fprintf ( output_unit, '  %d', element_node(order,element) - 1 );
    end
    fprintf ( output_unit, '\n' );
  end
%
%  VTK has a cell type 22 for quadratic triangles.  However, we
%  are going to strip the data down to linear triangles for now,
%  which is cell type 5.
%
  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'CELL_TYPES %d\n', element_num );

  if ( element_order == 4 )
    for element = 1 : element_num
        fprintf ( output_unit, '10\n' );
     
    end
%   elseif ( element_order == 6 )
%     for element = 1 : element_num
%       fprintf ( output_unit, '22\n' );
%     end
  end

  fprintf ( output_unit, '\n' );
  fprintf ( output_unit, 'POINT_DATA %d\n', node_num );
  fprintf ( output_unit, 'SCALARS pressure float\n' );
  fprintf ( output_unit, 'LOOKUP_TABLE default\n' );
  for node = 1 : node_num
    fprintf ( output_unit, '  %f\n', data(node,1) );
  end
  
  % for vector input of data
  if size(data,2)==3
  fprintf (output_unit,'\n');
  fprintf(output_unit,'VECTORS velocity float\n');
  for node=1:node_num
      fprintf(output_unit,'%f %f %f\n',[data(node,2:3),0.0])
  end
  end
  return
end
