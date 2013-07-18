function [next_ref_lev,cell_refine]=ref_difference_rule(next_ref_lev,cell_refine,...
    edg2poly,connectivity)

% number of cells flagged for refinement due to error indicator
n_cell_to_ref_ori = length(cell_refine);
% cells flagged for refinement due to error indicator
buffer = cell_refine;

% loop over cells flagged for refinement
while ~isempty(buffer)
    
    % pick a cell flagged for ref
    iel=buffer(1);
    buffer(1)=[];
        
    % proposed_ref_level_iel 
    proposed_ref_level_iel = next_ref_lev(iel);
    
    % get connectivity data
    g = connectivity{iel};
    
    % find the current neighbors to iel
    list_edge_p = find( edg2poly(:,2) == iel );
    Km = edg2poly(list_edge_p,1);
    % remove iel from list
    ind = find(Km==iel); Km(ind)=[];
    % remove boundary connections
    ind = find(Km<0);
    ind = sort(ind,'descend');
    for ii=1:length(ind)
        Km(ind(ii))=[];
    end
    list_edge_m = find( edg2poly(:,1) == iel );
    Kp = edg2poly(list_edge_m,2);
    % remove iel from list
    ind = find(Kp==iel); Kp(ind)=[];
    ind = find(Kp<0);
    % remove boundary connections
    ind = sort(ind,'descend');
    for ii=1:length(ind)
        Kp(ind(ii))=[];
    end
    K_others=[Kp ; Km];
    
    % loop over the neighbors of iel
    for k=1:length(K_others)
        level_diff = abs( next_ref_lev(K_others(k)) - proposed_ref_level_iel );
        switch level_diff
            case{0,1}
                % do nothing
            case{2}
                % flag k for refinement
                next_ref_lev(K_others(k)) = next_ref_lev(K_others(k)) + 1;
                buffer(end+1) = K_others(k);
                cell_refine(end+1) = K_others(k);
            otherwise
                error('level_diff should not be >2');
        end        
    end % end loop over the neighbors of iel
    
end % end while


% sanity check
aux = unique(cell_refine);
if(length(aux) ~= length(cell_refine))
    error('unique(cell_refine)');
end

% number of cells flagged for refinement due to error indicator + diff rule
n_cell_to_ref = length(cell_refine);

fprintf('Number of cells flagged for refinement due to error indicator %d\n',n_cell_to_ref_ori);
fprintf('Number of cells flagged for refinement due to err.ind + rule  %d\n',n_cell_to_ref);


return
end
