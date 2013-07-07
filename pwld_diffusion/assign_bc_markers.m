function edg2poly = assign_bc_markers(n_edge,edg2poly,edg2vert,vert,Lx,Ly)

% assign bc markers: LRBT = -( 10 20 30 40 )
for ied=1:n_edge
    % get K+ for that edge
    Kp = edg2poly(ied,2);
    % we want to loop only on BOUNDARY edges
    if(Kp>0), continue; end
    % get the 2 vertices associated with that edge
    P=edg2vert(ied,1:2);
    v=vert(P,:);
    x1=v(1,1); y1=v(1,2);
    x2=v(2,1); y2=v(2,2);
    % assign BC markers
    if(abs(x1)<1e-14 && abs(x2)<1e-14),
        edg2poly(ied,2)=-10; % left
    end
    if(abs(x1-Lx)<1e-14 && abs(x2-Lx)<1e-14),
        edg2poly(ied,2)=-20; % right
    end
    if(abs(y1)<1e-14 && abs(y2)<1e-14),
        edg2poly(ied,2)=-30; % bottom
    end
    if(abs(y1-Ly)<1e-14 && abs(y2-Ly)<1e-14),
        edg2poly(ied,2)=-40; % top
    end
end

return
end