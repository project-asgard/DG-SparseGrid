function val = cmerge( tcond,  true_value, false_value )
% val = cmerge( tcond,  true_value, false_value )
% simulate conditional assignment
% (tcond) ? (true_value) : (false_value)
% --------------------
if (tcond),
   val = true_value;
else
   val = false_value;
end;

