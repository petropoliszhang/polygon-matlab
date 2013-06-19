function rectangular_pwld()
clc; close all;
% clear all; close all; clc
%
% data
%
Lx=1; D=1; S=0; Q=0; Ly=1;
%
% numerical parameters
%
nx=3; ny=nx;
x=linspace(0,Lx,nx+1); y=linspace(0,Ly,ny+1);
nel=nx*ny;
ndof = 4*nel;
% 4---3   vertex anti-clockwise ordering,
% |   |
% 1---2
connectivity=zeros(nel,4);
for iel=1:nel
    skip = 4*(iel-1);
    i1 = skip + 1;
    i2 = skip + 2;
    i3 = skip + 3;
    i4 = skip + 4;
    connectivity(iel,:)=[i1 i2 i3 i4];
end
% DG vertex coordinates (they are duplicated for simplicity)
ind=0;
vert=zeros(ndof,2);
for iel=1:nel
    j = floor((iel-1)/nx) + 1;
    i = iel - (j-1)*nx;
%     [iel i j]
    vert(ind+1,1:2)=[x(i)   y(j)  ];
    vert(ind+2,1:2)=[x(i+1) y(j)  ];
    vert(ind+3,1:2)=[x(i+1) y(j+1)];
    vert(ind+4,1:2)=[x(i)   y(j+1)];
    ind = ind + 4;
end
% single coordinates
ind=0;
vert_grid=zeros((nx+1)*(ny+1),2);
for j=1:ny+1
    for i=1:nx+1
        ind=ind+1;
        vert_grid(ind,1:2)=[x(i) y(j)];
    end
end
% find relationship between vert and vert_grid
vert_link=zeros(ndof,1);
for i=1:ndof
    v=vert(i,1:2);
    for k=1:(nx+1)*(ny+1)
        if(norm( vert_grid(k,:)-v )<1e-12)
            if(vert_link(i)~=0), error('vert_link(i) should be 0'); end
            vert_link(i)=k;
        end
    end
end
% edge data
n_edge=0;
edg2poly=zeros(0,2);
edg_normal=zeros(0,2);
edg2vert=zeros(0,4);
for iel=1:nel
    elem=connectivity(iel,:);
    nedg=length(elem);
    elem(end+1)=elem(1);
    for i=1:nedg
        ed=elem(i:i+1);
        [edg2poly,edg2vert,n_edge]=is_edge_already_recorded(...
            ed,edg2poly,edg2vert,iel,vert_link,n_edge);
    end
end
clear vert_link; % not needed any longer
% compute edge normals
for ied=1:n_edge
    v1=vert(edg2vert(ied,1),:);
    v2=vert(edg2vert(ied,2),:);
    vec=v2-v1;
    vec=vec/norm(vec);
    edg_normal(ied,1:2)=[vec(2) -vec(1)];
end
% DG assemble volumetric terms
A = spalloc(ndof,ndof,9); b=zeros(ndof,1);
for iel=1:nel
    g=connectivity(iel,:);
    v=vert(g,:);
    [M,K,f,grad{iel}]=build_pwld_local_matrices(g,v);
    A(g(:),g(:)) = A(g(:),g(:)) + D*K +S*M;
    b(g(:)) = b(g(:)) + Q*f;
end
spy(A)
% DG assemble edge terms
%              
%           v2 ^  w1
%              |
%              |  n_ed  
%  Minus side  | --->    Plus side
%              |
%              |
%           v1 .  w2
%           
for ied=1:n_edge
    % get K-,K+ and their connectivities
    Kp = edg2poly(ied,2);
    Km = edg2poly(ied,1);
    % we want to loop only on INTERIOR edges
    if(Kp<0 | Km<0), continue; end
    % get the polygons' connectivities
    gp = connectivity(Kp,:);
    gm = connectivity(Km,:);
    % nbr of vertices in each poly
    nvp = length(gp);
    nvm = length(gm);
    % get normal
    ne = edg_normal(ied,1:2);
    % get edge vertices of K+ side
    W = edg2vert(ied,3:4);
    % get edge vertices of K- side
    V = edg2vert(ied,1:2);
    % get the local ID of the edge's first vertex in K+
    indp = find( gp == W(1) );
    if(length(indp) ~= 1), error('length(IDp) ~= 1'); end
    if(indp ~= nvp)
        IDp = [indp (indp+1) ];
    else
        IDp = [indp 1 ];
    end        
    % get the local ID of the edge's first vertex in K-
    indm = find( gm == V(1) );
    if(length(indm) ~= 1), error('length(IDm) ~= 1'); end
    if(indm ~= nvm)
        IDm = [indm (indm+1) ];
    else
        IDm = [indm 1 ];
    end
    % skipping indices
    skip_p = [(indp:nvp) (1:indp-1)];
    skip_m = [(indm:nvm) (1:indm-1)];

    % half-length current edge
    L=norm( diff(vert(V,:)) )/2;
    
    % build the local edge gradient matrices
    % [[phi]],{{D.grad(b).ne}} 
    %      = (phi+ - phi-)(D+ grad(b+).ne + D- grad(b-).ne)/2
    %      =  phi+ D+ grad(b+).ne/2 
    %       + phi+ D- grad(b-).ne/2
    %       - phi- D+ grad(b+).ne/2
    %       - phi- D- grad(b-).ne/2
    % compute row vector: n' * G (my_n is already retrieved as a 1x2 row
    % vector)
    % -/-
    rv_mm = ne * grad{Km}(:,:,indm);
    % edge matrix for this side
    cv=zeros(nv,1); cv(list_vert(1:2))=1;
    edgmat_mm = cv * rv_mm
    

    
end

spy(A)
% apply bc
bcnodes=1:nx+1;
bcval(1:length(bcnodes))=1;
bcnodes=[bcnodes (ny*(nx+1)+1:ndof)];
bcval(length(bcval)+1:length(bcnodes))=0;
% bcnodes=[bcnodes  (nx+2:nx+1:(ny-1)*(nx+1)+1)    ];
% bcnodes=[bcnodes ((nx+2:nx+1:(ny-1)*(nx+1)+1)+nx)];
for i=1:length(bcnodes)
    bd=bcnodes(i);
    A(bd,:)=0;
    b=b-A(:,bd)*bcval(i);
    A(:,bd)=0;
    A(bd,bd)=1;
    b(bd)=bcval(i);
end

%solve
z=A\b;

% plot
for iel=1:nel
    g=connectivity(iel,:);
    patch(vert(g,1),vert(g,2),z(g),z(g),'FaceColor','interp'); %,'LineStyle','none');
end
figure(2)
% plot on finer mesh
% 4---3   vertex anti-clockwise ordering,
% | c |
% 1---2
for iel=1:nel
    g=connectivity(iel,:);
    v=vert(g,:);
    c=mean(v);
    zc=mean(z(g));
    % alpha coef
    nv=length(g);
    alpha=1/nv;
    for i=1:nv
        i2=i+1; if(i==nv), i2=1; end
        xx=[ vert(g([i i2]),1); c(1)];
        yy=[ vert(g([i i2]),2); c(2)];
        zz=[ z(g([i i2])); zc];
        patch(xx,yy,zz,zz,'LineStyle','none');
    end
end

return
end
