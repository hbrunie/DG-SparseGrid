function Y = kron_mult2(A1,A2,X)
options.format = 'h'; options.round = 1; options.subnormal = 1;
chop([],options)
% Y = kron_mult2(A1,A2,X)
%
use_kron_multd = 0;
if (use_kron_multd),
 Acell{1} = chop(A1);
 Acell{2} = chop(A2);
 nkron = 2;
 Y = kron_multd(nkron,Acell,chop(X));
else
  Y = kronmult2(chop(A1),chop(A2),chop(X));
end;
