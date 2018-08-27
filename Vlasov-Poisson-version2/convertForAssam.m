function [stat] = convertForAssam()

files = what();

fNames = files.mat;

nF = numel(fNames);

threshold = 1e-12;

for f=1:nF
    
    varsToClear = {'A','Mat'};
    clear(varsToClear{:});
    
    load([fNames{f}]);
    
    if exist('Mat','var')
        A = Mat;
    end
     
    fTxt = [fNames{f},'.txt'];
    
    disp(fTxt);
    
    [m,n] = size(A);
    disp(['m x n : ', num2str(m),' x ', num2str(n)]);
    
    fileID = fopen(fTxt,'w');
    
    %A(find(abs(A)<threshold))=0;
    
    [iAll,jAll,vAll] = find(A);
    
    data = [iAll(:)';jAll(:)';vAll(:)'];
    
    [M,N] = size(data);
    disp(['N (number non-zeros): ',num2str( N)]);
    
    if N<5e5
    RC = rcond(full(A));
    disp(['rcond : ', num2str(RC)] );
    end
    
    NZF = nnz(A)/numel(A)*100;
    disp(['nnz/numel(A)*100 (non-zero %): ',num2str(NZF), '%'] );
    
    disp('  ');
    
    fprintf(fileID, "%10d \n", n);
    fprintf(fileID, "%10d \n", N);
%     fprintf(fileID, "rcond : %25.8e \n", RC );
%     fprintf(fileID, "nnz/numel(A)*100 : 25.8e \n", NZF );
    fprintf(fileID, "%10d %10d %25.8e \n", data);
    
    fclose(fileID);
    
end

end