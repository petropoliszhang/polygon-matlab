function edg_normal = compute_edge_normals(n_edge,edg2vert,vert)

% init array
edg_normal = zeros(n_edge,2);

% compute edge normals
for ied=1:n_edge
    v1=vert(edg2vert(ied,1),:);
    v2=vert(edg2vert(ied,2),:);
    vec=v2-v1;
    if(norm(vec)<1e-10)
        warning('norm(vec) too small'); 
    end
    vec=vec/norm(vec);
    edg_normal(ied,1:2)=[vec(2) -vec(1)];
end

return
end
