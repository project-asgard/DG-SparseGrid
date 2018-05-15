function [] = write_fval_to_file(fval,Lev,Deg,timestep)

LevDegStr = sprintf('lev-%i_deg-%i',Lev,Deg);
TimeStepStr = sprintf('-%4.4i',timestep);

if ~exist(['data'],'dir'),
   [status, msg, msgid] = mkdir('data');
   isok = (status == 1);
   if (~isok),
     error(sprintf('problem in creating directory %s,msg=%g,msgid=%g', ...
             'data', msg, msgid ));
   end;
end;

if ~exist(['data/' LevDegStr],'dir'), 
  [status, msg, msgid] = mkdir('data/',LevDegStr); 
  isok = (status == 1);
  if (~isok),
    error(sprintf('problem in creating director %s,msg=%g,msgid=%g', ...
            ['data/' LevDegStr], msg, msgid));
  end;
end;

options = '';
be_compatible_with_matlab = 0;
if (be_compatible_with_matlab),
  options = '-ascii';
end;

save(['data/' LevDegStr '/fval' TimeStepStr '.mat'],'fval', options);

end
