x_node = 0:0.01:1;
% [XX,YY,ZZ] = meshgrid(x_node);
[YY,ZZ,XX] = meshgrid(x_node);
EEval = pde.E(XX(:),YY(:),ZZ(:),time);
E1 = EEval(:,:,:,1);
E2 = EEval(:,:,:,2); 
E3 = EEval(:,:,:,3);
BBval = pde.B(XX(:),YY(:),ZZ(:),time);
B1 = BBval(:,:,:,1);
B2 = BBval(:,:,:,2); 
B3 = BBval(:,:,:,3);
vtkwrite('Max1.vtk','structured_grid',XX,YY,ZZ,...
    'vectors','E',E1,E2,E3,'vectors','B',B1,B2,B3);

pde = Maxwell2;
time = 1/3+1/4+1/4+1/4;
EEval = pde.E(XX(:),YY(:),ZZ(:),time);
E1 = EEval(:,:,:,1);
E2 = EEval(:,:,:,2); 
E3 = EEval(:,:,:,3);
BBval = pde.B(XX(:),YY(:),ZZ(:),time);
B1 = BBval(:,:,:,1);
B2 = BBval(:,:,:,2); 
B3 = BBval(:,:,:,3);
vtkwrite('Max2_t1over3-3.vtk','structured_grid',XX,YY,ZZ,...
    'vectors','E',E1,E2,E3,'vectors','B',B1,B2,B3);