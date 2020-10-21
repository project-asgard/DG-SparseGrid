function f = real_space_solution_at_coordinates_irregular(pde,wavelet_coeffs,coordinates)

for p=1:numel(coordinates(:,1))
    
    f(p) = real_space_solution_at_coordinates(pde,wavelet_coeffs,coordinates(p,:));
    
end

figure
x = coordinates(:,1);
y = coordinates(:,2);
tri = delaunay(x,y);
trisurf(tri,x,y,f);

end
