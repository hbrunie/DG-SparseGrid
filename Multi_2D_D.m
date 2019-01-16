function [f_rSpace]=Multi_2D_D(Meval_D,f_wSpace,HASHInv,pde)

%%
% Apply dimension many inverse wavelet transforms to the solutionn vector
% f_wSpace -> f_rSpace (wavelet space -> real space)
% via Kron(tmp1,tmp2,...,tmpD)*f
% Note: the output grid of this transform is specified in matrix_plot_D.m

dimensions = pde.dimensions;

nDims = numel(dimensions);
nHash = numel(HASHInv);

deg = dimensions{1}.deg; % TODO : generalize to deg_D
lev = dimensions{1}.lev; % TODO : generalize to lev_D

%%
% Catch for TODO

if nDims>1
    for d=2:nDims
        assert(dimensions{d}.lev==dimensions{d}.lev);
        assert(dimensions{d}.deg==dimensions{d}.deg);
    end
end

%%
% TODO : surely this size depends on the parameters in the matrix_plot_D
% routine? i.e., the grid upon which we decide to evaluate the solution?
dof_1D_FG = deg*2^(lev); 

%%
% The real space vector is on the full-grid size? Where is this decided?
% (matrix_plot_D.m I think)
fnew = sparse(dof_1D_FG^nDims,1);
f_rSpace = sparse(dof_1D_FG^nDims,1);

%% 
% To be removed
if nDims==2
    A = Meval_D{1};
    B = Meval_D{2};
    f = f_wSpace;
end

for i=1:nHash
    
    ll=HASHInv{i};
    
    %%
    % Retain 2D version for checking but don't use
    if nDims==2
        I1=HASHInv{i}(5);
        I2=HASHInv{i}(6);
        
        index_I1=[(I1-1)*deg+1:I1*deg];
        index_I2=[(I2-1)*deg+1:I2*deg];
        
        Index = deg^2*(i-1)+1:deg^2*i;
        
        tmp=kron(...
            A(:,index_I1),...
            B(:,index_I2) ...
            )*f(Index(:));
        fnew=fnew+tmp;
    end
    
    %%
    % kron(tmp1,tmp2,...,tmpD)*fwav
    
    clear kronMatList;
    for d=1:nDims
        ID = ll(nDims*2+d); % TODO : Check if this indexing correct for D != 2?
        index_D = [(ID-1)*deg+1 : ID*deg];
       
        %%
        % Assert match the 2D case
        if nDims==2
            if d==1
                assert(norm(index_D-index_I1)==0);
            end
            if d==2
                assert(norm(index_D-index_I2)==0);
            end
        end
        
        thisMat = Meval_D{d};
        thisMat1 = thisMat(:,index_D);
        
        %%
        % Assert match the 2D case
        if nDims==2
            if d==1
                assert(norm(thisMat1-A(:,index_I1))==0);
            end
            if d==2
                assert(norm(thisMat1-B(:,index_I2))==0);
            end
        end
        
        kronMatList{d} = thisMat1; % Build kron list
    end
   
    element_ii = deg^nDims*(i-1)+1:deg^nDims*i;
    
    X = f_wSpace(element_ii);
    Y = kron_multd(nDims,kronMatList,X);
    
    tol = 1e-15;
    if nDims==2
        assert(norm(Y-tmp)<tol);
    end
    
    f_rSpace = f_rSpace + Y;
    
end

if nDims==2
    assert(norm(fnew-f_rSpace)<tol);
end

end



