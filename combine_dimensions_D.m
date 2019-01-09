function [fval] = combine_dimensions_D(fD,ft,HASHInv,pde)

% Combine (via kron product) a set of 1D multiwavelet transforms to form
% the higher D sparse-grid multiwavelet representation.

% ft is the time multiplier.

nDims = numel(fD);
nHash = numel(HASHInv);

Deg = pde.dimensions{1}.deg; % TODO Need to generalize this deg_D

fval = sparse(Deg^nDims * nHash,1);

for i=1:nHash
    
    ll=HASHInv{i};
    
    %%
    % Kron product approach
    
    A = 1;
    for d=1:nDims
        ID = ll(nDims*2+d); % TODO : Check if this indexing correct for D != 2?
        index_D = [(ID-1)*Deg+1 : ID*Deg];
        f = fD{d};
        ftmp = f(index_D);
        A = kron(A,ftmp);
    end
    
    B = ft;
    
    tmp = A * B;
    
    Index = Deg^nDims*(i-1)+1:Deg^nDims*i;
    fval(Index,1) = fval(Index,1) + tmp(:);
    
end

end


