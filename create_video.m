function create_video()

% List output files

folder = './';
s = what(folder);

list = s.mat;

N = numel(list);

% Initalize video file

v = VideoWriter('fk6.avi');
open(v);

for n=1:N
    
    load(list{n});
   
    contourf(f2d,'EdgeColor','none');
    %caxis([range2n +range2]);
    title('f');
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end

close(v);

end