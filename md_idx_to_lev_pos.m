function [lev_vec, pos_vec] = md_idx_to_lev_pos (num_dims, max_lev, idx)

% disp(['md_idx_to_lev_pos max_lev : ', num2str(max_lev)]);

stride = 2^max_lev;

sd_idx = zeros(num_dims,1);

for d=1:num_dims
    
    sd_idx(d) = fix(mod((idx-1)/stride^(d-1),stride));
    
    [lev,pos] = sd_idx_to_lev_pos(sd_idx(d)+1);
    lev_vec(d) = lev;
    pos_vec(d) = pos;
    
end

assert(min(lev_vec)>=0);
assert(min(pos_vec)>=0);

%assert(max(lev_vec)<=max_lev,'not sure why this doesnt work, but just increase max_lev');

end


% function [lev_vec, pos_vec] = md_idx_to_lev_pos (num_dims, max_lev, idx)
% 
% if num_dims == 1
%     
%     [lev,pos] = sd_idx_to_lev_pos(idx);
%     
%     lev_vec(1) = lev;
%     pos_vec(1) = pos;
%     
% elseif num_dims == 2
%     
%     stride = prod(ones(num_dims-1,1)*(2^max_lev))
%     
%     sd_idx_1 = mod(idx-1,stride);
%     sd_idx_2 = (idx-1)/stride;
%     
%     [lev1,pos1] = sd_idx_to_lev_pos(sd_idx_1+1);
%     [lev2,pos2] = sd_idx_to_lev_pos(sd_idx_2+1);
%     
%     lev_vec(1) = lev1;
%     lev_vec(2) = lev2;
%     
%     pos_vec(1) = pos1;
%     pos_vec(2) = pos2;
%     
% elseif num_dims == 3
%     
%     stride = 2^max_lev;
%     
%     sd_idx = zeros(num_dims);
%     for d=1:num_dims
%         
%         sd_idx(1) = mod(idx-1,stride);
%         sd_idx(2) = mod((idx-1)/stride,stride);
%         sd_idx(3) = mod((idx-1)/stride^2,stride);
%         
%     end
%     
%     [lev1,pos1] = sd_idx_to_lev_pos(sd_idx_1+1);
%     [lev2,pos2] = sd_idx_to_lev_pos(sd_idx_2+1);
%     [lev3,pos3] = sd_idx_to_lev_pos(sd_idx_3+1);
%     
%     lev_vec(1) = lev1;
%     lev_vec(2) = lev2;
%     lev_vec(3) = lev3;
%     
%     pos_vec(1) = pos1;
%     pos_vec(2) = pos2;
%     pos_vec(3) = pos3;
%     
% end
% 
% end
