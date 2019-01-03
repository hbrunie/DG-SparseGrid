function A_encode=GlobalMatrixSG(A,B,HASH,Lev,Deg)
% Global Matrix construction from Matrices A and B
% by Hash Table
% This code is to use the new approach for generating global matrix
% only works for
% compression = 3
% By Lin, 11/14/2018

global hash_format

Dim = 2;
ichoice = 1;


ComLev = perm_leq( Dim, Lev);
ComLevIndex = [];

for i = 1:size(ComLev,1)
    
    Lev_loc = ComLev(i,:);
    
    for d = 1:Dim
        Cell_loc{d}=[0:2^max(Lev_loc(d)-1,0)-1];
        index = LevCell2index(Lev_loc(d),Cell_loc{d});
        ComLevIndex{i}.lid{d} = sort(index);
    end
    
    if ichoice == 1
        [Cel1,Cel2] = meshgrid(Cell_loc{:}); % need more work to generalize d
                Cel1 = Cel1';
                Cel2 = Cel2';
    else
        nz1=size(Cell_loc{1},2);
        nz2=size(Cell_loc{2},2);
        Cel1 = zeros(nz1*nz2,1);
        Cel2 = zeros(nz1*nz2,1);
        for k = 1:nz1
            Cel1(k:nz2:end)=Cell_loc{1}(k);
        end
        for k = 1:nz2
            Cel2((k-1)*nz1+1:k*nz1) = Cell_loc{2}(k);
        end
    end
    
    key = zeros(numel(Cel1),2*Dim);
    key(1:end,1:Dim) = repmat(Lev_loc,numel(Cel1),1);
    key(:,Dim+1:2*Dim) = [Cel1(:),Cel2(:)];
    
    %     i
    index = zeros(size(key,1),1);
    % find the index from HashTable
    for j =1:size(key,1)
        index(j) = HASH.(sprintf(hash_format,key(j,:)));
    end
    
    ComLevIndex{i}.gid = index;
    clear key Cel1 Cel2 index
end
count = 1;

for i = 1:size(ComLev,1)
    for j = 1:size(ComLev,1)
        IndexI = ComLevIndex{i}.gid;
        sizeI = size(IndexI,2);
        IndexJ = ComLevIndex{j}.gid;
        sizeJ = size(IndexJ,2);
        
        for d = 1:Dim
            index_I{d} = ComLevIndex{i}.lid{d};
            index_J{d} = ComLevIndex{j}.lid{d};
        end
        
        tmpA=A( Deg*(index_I{1}(1)-1)+1:Deg*(index_I{1}(end)),...
            Deg*(index_J{1}(1)-1)+1:Deg*(index_J{1}(end)) );
        rA = size(tmpA,1);
        cA = size(tmpA,2);
        
        tmpB=B( Deg*(index_I{2}(1)-1)+1:Deg*(index_I{2}(end)),...
            Deg*(index_J{2}(1)-1)+1:Deg*(index_J{2}(end)) );
        rB = size(tmpB,1);
        cB = size(tmpB,2);        
        
        A_encode{count}.A1=tmpA;
        A_encode{count}.A2=tmpB;
        
%         if size(tmpB,1)/Deg>1 || size(tmpB,2)/Deg>1 ||
%         if i == 10 || j==10
%             1111
%         end
        
        
        tmp = Deg^2*(IndexI(:)-1)+[1:Deg^2];
        tmp = tmp';
        rIndex = kron_split( rA,  Deg*ones(rB/Deg,1) );
%         rIndex = rIndex + Deg^2*[0:rA/Deg-1];
        
        IndexI = tmp(:);
        IndexI(rIndex(:)) = IndexI;
        
        
        tmp = Deg^2*(IndexJ(:)-1)+[1:Deg^2];
        tmp = tmp';
        cIndex = kron_split( cA,  Deg*ones(cB/Deg,1) );
%         cIndex = cIndex + Deg^2*[0:cA/Deg-1];
        IndexJ = tmp(:);
        IndexJ(cIndex(:)) = IndexJ;
        
        A_encode{count}.IndexI = IndexI;
        A_encode{count}.IndexJ = IndexJ;
        
        
        count = count+1;
    end
end

