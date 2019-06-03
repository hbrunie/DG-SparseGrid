function [elements,elementsIDX] = element_table(pde,opts)

num_dimensions = numel(pde.dimensions);
is_sparse_grid = strcmp( opts.gridType, 'SG');

%%
% Setup element table as a collection of sparse vectors to
% store the lev and cell info for each dim. 

%%
% set the maximum number of elements we can refine to

num_elements_max = 1;
lev  = 0;

for d=1:num_dimensions
    
%     this_num_elements = uint64(2)^pde.dimension{d}.lev;
    this_num_elements = uint64(2)^pde.maxLev;

    num_elements_max  = num_elements_max * this_num_elements;
    
    if pde.dimensions{d}.lev > lev
        lev = pde.dimensions{d}.lev;
    end
    
end
num_elements_max = double(num_elements_max); % This number is HUGE

%%
% allocate the sparse element table members

elements.lev_p1     = sparse (num_elements_max, num_dimensions); % _p1 is for "plus 1" sinse sparse cannot accpet 0
elements.pos_p1     = sparse (num_elements_max, num_dimensions);
elements.node_type  = sparse (num_elements_max, 1);

%%
% Get the combinations of levels

if (is_sparse_grid)
   ptable = perm_leq( num_dimensions, lev );
else
   ptable = perm_max( num_dimensions, lev );
end

%%
% Remove lev values not allowed due to rectangularity (yes, it is a word)

keep = ones(size(ptable,1),1);
for i=1:size (ptable, 1)
    for d=1:num_dimensions
        if (ptable(i,d) > pde.dimensions{d}.lev)
            keep(i) = 0;
        end
    end
end

ptable = ptable(find(keep),:);

%%
% compute the number of cells for each combination of levels

ncase = size(ptable,1);
isize = zeros(1,ncase);

levels = zeros(1,num_dimensions);
ipow = zeros(1,num_dimensions);

for icase=1:ncase
   levels(1:num_dimensions) = ptable(icase,1:num_dimensions);
   ipow(1:num_dimensions) = 2.^max(0,levels(1:num_dimensions)-1);
   isize(icase) = prod( max(1,ipow) );
end

% ---------------------------
% prefix sum or cumulate sum
% note   matlab
%   cumsum( [2 3 4 ]) 
%   returns
%           [2 5 9 ]
%
%  a a  b b b  c c c c
%  istart contains
%  1    3      6       10
% ---------------------------
istartv = zeros(1,ncase+1);
istartv(1) = 1;
istartv(2:(ncase+1)) = cumsum(isize);
istartv(2:(ncase+1)) = istartv(2:(ncase+1)) + 1;

for icase=1:ncase
  istart = istartv(icase);
  iend   = istartv(icase+1)-1;

  levels(1:num_dimensions) = ptable(icase,1:num_dimensions);
  index_set = lev_cell_to_singleD_index_set( levels(1:num_dimensions) );
  
  for i=istart:iend
      
     icells = index_set(i-istart+1,:);
     
     %%
     % Store the index into the element table for this element
     
     element_idx = lev_cell_to_element_index(pde,levels,icells);
     elementsIDX(i) = element_idx;
     
     %%
     % Set the lev and cell coordinate for each dim
     
     elements.lev_p1(element_idx,:) = levels+1; % NOTE : have to start lev  index from 1 for sparse storage
     elements.pos_p1(element_idx,:) = icells+1; % NOTE : have to start lev  index from 1 for sparse storage
     
     %%
     % Set the element type ("leaf" or "internal") by checking if the sum
     % over level across dimensions equals the sparse grid selection limit
     % from perm_leq (pde, lev)
     
     elements.node_type(element_idx) = 1; % 'internal'; % Internale nodes will not be checked for refinement.
     if sum (elements.lev_p1(element_idx,:)-1) == lev
         elements.node_type(element_idx) = 2; % 'leaf'; % Leaf nodes are checked for refinement
     end

  end
end


end
