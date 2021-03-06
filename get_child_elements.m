function [child_elements_lev_vec, child_elements_pos_vec, cnt] = ...
    get_child_elements (lev_vec, pos_vec, max_lev, refinement_method)


%%
% Takes lev and pos vector as input, returns the same for a list of
% daughter elements which also obey one of the sparse adaption rules.

num_dimensions = numel(lev_vec);

assert (num_dimensions == numel (pos_vec) );
assert (min(lev_vec)>=0);
assert (min(pos_vec)>=0);

debug = 0;

child_elements_lev_vec = [];
child_elements_pos_vec = [];

%%
% Available refinement methods
%
% In this 2D example, assume we want to refine the (2,b) element and do so
% in all directions (dimensions)
%
%         1                a
%        / \              / \
%       2   3            b   c
%      / \              / \
%     4   5            d   e
%
% Method 1 (david)
% 
% number of new elements = 2 * num_dimension
%
% (4,b)
% (5,b)
% (2,d)
% (2,e)
%
% Method 2 (david + lin)
%
% number of new elements = 2 * num_dimension + 2^num_dimensions
%
% (4,b)
% (5,b)
% (2,d)
% (2,e)
% (4,d)
% (4,e)
% (5,d)
% (5,e)

if ~exist('refinement_method','var') || isempty(refinement_method)
    refinement_method = 1;
end

cnt = 0;
if refinement_method == 1 || refinement_method == 2
    
    for d=1:num_dimensions                   
                    
        %%
        % First daughter
        
        new_elem_lev_vec = lev_vec;
        new_elem_pos_vec = pos_vec;
                
        if new_elem_lev_vec(d)+1 < max_lev
            
            new_elem_lev_vec(d) = new_elem_lev_vec(d)+1;
            new_elem_pos_vec(d) = new_elem_pos_vec(d)*2; % Assumes pos starts at 0
            
            child_elements_lev_vec(cnt+1,:) = new_elem_lev_vec;
            child_elements_pos_vec(cnt+1,:) = new_elem_pos_vec; % Assumes pos starts at 0
            
            if child_elements_lev_vec(cnt+1,d) < 0
                disp('l');
            end
            assert(child_elements_pos_vec(cnt+1,d) >= 0);
            assert(child_elements_lev_vec(cnt+1,d) >= 0);
            
            cnt = cnt + 1;
            
            new_lev_1D{d} = new_elem_lev_vec(d);
            new_pos_1D{d} = new_elem_pos_vec(d);
            
        end
        
        
        %%
        % If level is 1 or deeper add a second daughter
        
        if lev_vec(d) >= 1
                        
            new_elem_lev_vec = lev_vec;
            new_elem_pos_vec = pos_vec;
            
            if new_elem_lev_vec(d)+1 < max_lev
                
                new_elem_lev_vec(d) = new_elem_lev_vec(d)+1;
                new_elem_pos_vec(d) = new_elem_pos_vec(d)*2+1; % Assumes pos starts at 0
                
                child_elements_lev_vec(cnt+1,:) = new_elem_lev_vec;
                child_elements_pos_vec(cnt+1,:) = new_elem_pos_vec; % Assumes pos starts at 0
                
                assert(child_elements_pos_vec(cnt+1,d) >= 0);
                assert(child_elements_lev_vec(cnt+1,d) >= 0);
                
                cnt = cnt + 1;
                
                new_lev_1D{d} = [new_lev_1D{d},new_elem_lev_vec(d)];
                new_pos_1D{d} = [new_pos_1D{d},new_elem_pos_vec(d)];
                
            end
            
        end
        
    end
    
end

if refinement_method == 2
    
    if num_dimensions == 1
        
        %%
        % Do nothing as method 2 == method 1 for 1D
        
    end
    
    if num_dimensions == 2
        
        for i=1:numel(new_lev_1D{1})
            for j=1:numel(new_lev_1D{2})
                
                tmp_lev(1) = new_lev_1D{1}(i);
                tmp_lev(2) = new_lev_1D{2}(j);
                
                tmp_pos(1) = new_pos_1D{1}(i);
                tmp_pos(2) = new_pos_1D{2}(j);
                
                if max(tmp_lev) < max_lev
                    
                    child_elements_lev_vec(cnt+1,1:2) = tmp_lev;
                    child_elements_pos_vec(cnt+1,1:2) = tmp_pos;
                    
                end
                
                cnt = cnt + 1;
                
            end
        end
        
    end
    
    if num_dimensions == 3
        
        for i=1:numel(new_lev_1D{1})
            for j=1:numel(new_lev_1D{2})
                for k=1:numel(new_lev_1D{3})
                    
                    
                    tmp_lev(1) = new_lev_1D{1}(i);
                    tmp_lev(2) = new_lev_1D{2}(j);
                    tmp_lev(3) = new_lev_1D{3}(k);
                                   
                    tmp_pos(1) = new_pos_1D{1}(i);
                    tmp_pos(2) = new_pos_1D{2}(j);
                    tmp_pos(3) = new_pos_1D{3}(k);
                    
                    if max(tmp_lev) < max_lev
                        child_elements_lev_vec(cnt+1,1:3) = tmp_lev;
                        child_elements_pos_vec(cnt+1,1:3) = tmp_pos;
                    end
                    
                    cnt = cnt + 1;
                    
                end
            end
        end
        
    end
    
    if num_dimensions > 3
        
        %%
        % TODO
        
        assert(1==2);
        
    end
    
end

if cnt>0
    assert(min(child_elements_lev_vec(:))>=0);
    assert(min(child_elements_pos_vec(:))>=0);
end

end
