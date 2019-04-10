function [total_flops, max_flops, min_flops] = flops_batch_list( batch_list )
% [total_flops,max_flops, min_flops] = flops_batch_list( batch_list )
%
% estimate total  amount of work
%
total_flops = 0;
max_flops = 0;
min_flops = 10^18;

nbatch = batch_list.nbatch;
for ibatch=1:nbatch,
  mm = batch_list.mlist(ibatch);
  nn = batch_list.nlist(ibatch);
  kk = batch_list.klist(ibatch);

  flops = (2.0*mm)*nn*kk;
  max_flops = max( max_flops, flops);
  min_flops = min( min_flops, flops);
  total_flops = total_flops + flops;

end;

  

