function Y = kron_multd_full( nkron, Acell, X )

options.format = 'h'; options.round = 1; options.subnormal = 1;
chop([],options)
Y = 1;
for d=1:nkron
    
    Y = kron(Y,Acell{d});
    
end

Y = chop(chop(Y) * chop(X));

end
