function Y = krond(nkron, Acell)
% Y = krond(nkron, Acell)
% n-dimensional version of kron
%
% Y = kron( Acell{1}, Acell{2}, ..., Acell{nkron} )
%
options.format = 'h'; options.round = 1; options.subnormal = 1;
chop([],options)
if (nkron == 1),
    Y = chop(Acell{1});
    return;
end

if (nkron == 2),
    Y = chop(kron( chop(Acell{1}), chop(Acell{2}) ));
    return;
end

Y = chop(kron(krond( nkron-1, Acell), chop(Acell{nkron})));
end
