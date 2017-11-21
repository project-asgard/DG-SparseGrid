% 2D test
% clear

% n=2;k=2;
count=1;
database=struct();

for n1=0:n
    for i1=0:max(0,2^max(0,n1-1)-1)
        
        
        for n2=0:n-n1
            for i2=0:max(0,2^max(0,n2-1)-1)
                
                for k1=0:k-1
                    for k2=0:k-1
                        
                        key=[n1,n2,i1,i2,k1,k2];
                        database.(sprintf('i%g_',key))=count;
                        inv{count}=key;
                        
                        count=count+1;
                        
                        
                    end
                end
                
            end
        end
        
    end
end

dof_sparse=count-1;

% key=[0,0,0,0,0,0];
% database.(sprintf("i%g_",key))

