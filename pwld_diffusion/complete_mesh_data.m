function [n_edge,edg2poly,edg2vert,edg_perp ] = complete_mesh_data( nel,ndof,n_grid_vert,vert,vert_grid,connectivity )
% complete_mesh_data
global verbose

% find relationship between vert (DG) and vert_grid
vert_link=zeros(ndof,1);
% % % for kicks:
% % for i=n_grid_vert:-1:2
% %     vert_grid(i,:)=[];
% % end
% % n_grid_vert=1;

for i=1:ndof
    v=vert(i,1:2);
    % get difference
    aux = vert_grid - kron(ones(n_grid_vert,1),v);
    % get vector of norm
    aa=sqrt(aux(:,1).^2+aux(:,2).^2);
    ind=find(aa<1e-12);
    len=length(ind);
    if(len >1),
        v
        ind
        len
        vert_grid(ind,:)
        error('vert_link');
    end
    if(len==0),
        if(verbose)
            v
            ind
            len
            vert_grid(ind,:)
            warning('missing vert_link');
        end
        n_grid_vert = n_grid_vert + 1;
        vert_grid(n_grid_vert,:) = v;
        ind=n_grid_vert;
    end
    vert_link(i)=ind(1);
end

% edge data
n_edge=0;
edg2poly=zeros(0,2);
edg2vert=zeros(0,4);
for iel=1:nel
    g=connectivity{iel}(:);
    nedg=length(g);
    g(end+1)=g(1);
    for i=1:nedg
        ed=g(i:i+1);
        [edg2poly,edg2vert,n_edge]=is_edge_already_recorded(...
            ed,edg2poly,edg2vert,iel,vert_link,n_edge);
    end
end
clear vert_link; % not needed any longer

% check orientation, centroid
tot_area=0;
centroid_outside=0;
centroid_on_edge=0;
for iel=1:nel
    g=connectivity{iel}(:);
    xx=vert(g,1); yy=vert(g,2);
    
    % check orientation, verify area
    [or,ar] = polyorient(xx,yy);
    if(or~=1)
        iel
        g
        [xx yy]
        or
        ar
        error('orientation problem');
    end
    tot_area=tot_area+ar;
    
    % check ''centroid''
    xc=mean(xx); yc=mean(yy);
    [in on]=inpolygon(xc,yc,xx,yy);
    % find in points
    ind_in=find(in==1);
    % find if these 'in' points are on the poly edge
    in_and_on = on(ind_in);
    ind_in_and_on=find(in_and_on==1);
    if(~isempty(ind_in_and_on))
        if(verbose), warning('centroid is on poly edge'); end
        centroid_on_edge = centroid_on_edge + 1;
    end
    % find out points
    ind_out=find(in==0);
    if(~isempty(ind_out))
        if(verbose)
            iel
            g
            warning('centroid is OUTSIDE of poly');
        end
        centroid_outside = centroid_outside + 1;
    end
    
end
fprintf('total area read in geom = %g \n',tot_area);
fprintf('%d centroid(s) outside out of %d cells \n',centroid_outside,nel);
fprintf('fraction of centroids outside of polygon = %g \n',centroid_outside/nel);
fprintf('fraction of centroids on edge of polygon = %g \n',centroid_on_edge/nel);

% compute h_perp
edg_perp=zeros(n_edge,2); % j=1 for Km, j=2 for Kp
edg_perp_=zeros(n_edge,2); % j=1 for Km, j=2 for Kp
for ied=1:n_edge
    % get K-,K+
    Kp = edg2poly(ied,2);
    Km = edg2poly(ied,1);
    if(Km<=0 || Kp==0), error('Km<=0 or Kp==0'); end
    
    % get the Km polygon connectivity and vertices
    gm = connectivity{Km}(:);
    xx=vert(gm,1); yy=vert(gm,2);
    % get area
    [or,poly_area] = polyorient(xx,yy);
    % get perimeter
    xx=[xx; xx(1)]; yy=[yy; yy(1)];
    xx=diff(xx); yy=diff(yy);
    perim=sum(sqrt( xx.^2 + yy.^2 ));
    % get edge vertices of K- side
    V = edg2vert(ied,1:2);
    % length current edge
    Le = norm( diff(vert(V,:)) );
    % nbr of vertices in Km poly
    nvm = length(gm);
    if nvm <3,
        error('the # of vertices cannot be <3');
    elseif nvm==3
        edg_perp(ied,1) = 2*poly_area/Le;
    elseif nvm==4
        edg_perp(ied,1) = poly_area/Le;
    elseif  mod(nvm,2)==0
        edg_perp(ied,1) = 4*poly_area/perim;
    else
        edg_perp(ied,1) = 2*poly_area/perim + sqrt(2*poly_area/nvm/sin(2*pi/nvm));
    end
    % new option
    if(nvm > 3)
        X1 = vert(V(1),:);
        dx = diff(vert(V,:));
        hbot = 1e99;
        % get points of Km that are not on that edge
        g_diff = setdiff(gm,V);
        if(length(g_diff) ~= nvm-2), error('nvm-2 hbot'); end
        for iv=1:length(g_diff)
            X0 = vert(g_diff(iv),:);
            dist = abs( det([(X1-X0)' dx']) ) / Le;
            hbot = min(hbot,dist);
        end
        edg_perp_(ied,1) = hbot;
        if(hbot>edg_perp(ied,1))
%             edg_perp(ied,1) = (edg_perp(ied,1) + hbot)/2;
%             edg_perp(ied,1) = hbot;
        end
    end
    
    if(Kp>0)
        gp = connectivity{Kp}(:);
        xx=vert(gp,1); yy=vert(gp,2);
        % get area
        [or,poly_area] = polyorient(xx,yy);
        % get perimeter
        xx=[xx; xx(1)]; yy=[yy; yy(1)];
        xx=diff(xx); yy=diff(yy);
        perim=sum(sqrt( xx.^2 + yy.^2 ));
        % nbr of vertices in Kp poly
        nvp = length(gp);
        if nvp <3
            error('the # of vertices cannot be <3');
        elseif nvp==3
            edg_perp(ied,2) = 2*poly_area/Le;
        elseif nvp==4
            edg_perp(ied,2) = poly_area/Le;
        elseif  mod(nvp,2)==0
            edg_perp(ied,2) = 4*poly_area/perim;
        else
            edg_perp(ied,2) = 2*poly_area/perim + sqrt(2*poly_area/nvp/sin(2*pi/nvp));
        end
        % new option
        % get edge vertices of K+ side
        W = edg2vert(ied,3:4);
        if(nvp > 3)
            X1 = vert(W(1),:);
            dx = diff(vert(W,:));
            hbot = 1e99;
            % get points of Km that are not on that edge
            g_diff = setdiff(gp,W);
            if(length(g_diff) ~= nvp-2), error('nvp-2 hbot'); end
            for iv=1:length(g_diff)
                X0 = vert(g_diff(iv),:);
                dist = abs( det([(X1-X0)' dx']) ) / Le;
                hbot = min(hbot,dist);
            end
            edg_perp_(ied,2) = hbot;
            if(hbot>edg_perp(ied,2))
%                 edg_perp(ied,2) = (edg_perp(ied,2) + hbot)/2;
%                 edg_perp(ied,2) = hbot;
            end
        end
    end
    
end

return
end

