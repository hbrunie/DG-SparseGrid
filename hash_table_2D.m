function [forwardHash,inverseHash] = hash_table_2D(Lev,Dim,gridType)
% [forwardHash,inverseHash] = HashTable(Lev,Dim,gridType)
%-------------------------------------------------
% Generate 2D Hash Table s.t n1+n2<=Lev
% Input: Lev:: Level information
%        Dim:: Dimensionality
% Output: forwardHash:: HashTable
%         inverseHash:: Inverse Looking up for Hash
% Major Change:: ignoring the Deg from mesh
% Adding the 1D index into HashTable with
%   (Lev_1D,Cell_1D)->Index_1D
% so the inv = (Lev_1,Lev_2,Cell_1,Cell_2,Index_1,Index_2)
%        key = [Lev_1,Lev_2,Cell_1,Cell_2]
%-------------------------------------------------

if ~exist('gridType','var') || isempty(gridType)
    gridType = 'SG'; % Set default gridType to SG
end

% global hash_format
% 
% % Specifies the number of allowable integers in the elements of the hash key
% % If more are needed, i.e., > 99, then change to 'i%3.3i_'.
% 
% hash_format =  'i%04.4d';

%%
% Set the hash format in one place 
set_hash_format

count=1;
forwardHash = struct(); % Empty struct array
inverseHash = {}; % Empty cell array

for n1=0:Lev
    if strcmp(gridType,'FG')
        Lev_tmp = Lev;
    else
        Lev_tmp = Lev-n1;
    end
    for n2=0:Lev_tmp


            for i1=0:max(0,2^max(0,n1-1)-1)
            for i2=0:max(0,2^max(0,n2-1)-1)
                
                key=[n1,n2,i1,i2];
                forwardHash.(sprintf(hash_format,key)) = count;
                
                % Linearize the heirarchial multi-index for each dimension.
                
                index_1 = lev_cell_to_1D_index(n1,i1);
                index_2 = lev_cell_to_1D_index(n2,i2);
                
                inverseHash{count} = [key,index_1,index_2];
                
                count = count+1;
            end
            end
        
end
end

% Add some other useful information to the forwardHash struct

forwardHash.Lev = Lev;
forwardHash.Dim = Dim;

end
