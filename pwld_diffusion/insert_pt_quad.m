function [ old_connectivity,old_vert ] = insert_pt_quad( ed,K_others,old_connectivity,old_vert,cell_refine,new_pt )
% insert_pt_quad

insert_fn = @(val, x, pos)cat(2,  x(1:pos-1), val, x(pos:end));

Kbot=[];
for i=1:length(K_others)
    g_neigh = old_connectivity{K_others(i)};
%     ed
%     g_neigh
    combined = [ ed ; g_neigh ];
    if( length(unique(combined)) == 4 ) % quads only!!!
        Kbot = K_others(i);
        continue
    end
end

% there is a neighbor
if(~isempty(Kbot))
    ind = find( cell_refine == Kbot );
    if(~isempty(ind))
        % nothing to do, that neighbor cell is listed to be refined anyway
        if(length(ind)~=1), error(' Kbot duplicate?'); end
    else
        % that neighbor cell is NOT listed for refinement, add edge point
        %         i1 = find( g_neigh == ed(2) );
        i2 = find( g_neigh == ed(1) );
        n_old_vert = length(old_vert(:,1));
        g_neigh = insert_fn(n_old_vert+1, g_neigh, i2);
        old_connectivity{Kbot} = g_neigh;
        old_vert(n_old_vert+1,:) = new_pt;
    end
end

return
end

