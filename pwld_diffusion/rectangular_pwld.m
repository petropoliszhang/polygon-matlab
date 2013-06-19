function rectangular_pwld()
clc; close all;
% clear all; close all; clc
%
% data
%
Lx=1; c_diff=1; sigma_a=0; S_ext=0; Ly=1;
%
% numerical parameters
%
nx=3; ny=nx;
x=linspace(0,Lx,nx+1); y=linspace(0,Ly,ny+1);
nel=nx*ny;
i_mat=ones(nel,1);
ndof = 4*nel;
C_pen=2;
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
    mat = i_mat(iel);
    [M,K,f,grad{iel}]=build_pwld_local_matrices(g,v);
    A(g(:),g(:)) = A(g(:),g(:)) + c_diff(mat)*K +sigma_a(mat)*M;
    b(g(:)) = b(g(:)) + S_ext(mat)*f;
end
% spy(A)
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
m1d=[2 1 ; 1 2]/6;
m1d_mod=[1 2; 2 1]/6;
for ied=1:n_edge
    % get K-,K+ and their connectivities
    Kp = edg2poly(ied,2);
    Km = edg2poly(ied,1);
    % we want to loop only on INTERIOR edges
    if(Kp<0 || Km<0), continue; end
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
    skip_pr= [(indp:-1:1) (nvp:-1:indp+1)];
    skip_mr= [(indm:-1:1) (nvm:-1:indm+1)];

    % length current edge
    Le = norm( diff(vert(V,:)) );
    % material properties
    Dp = c_diff(i_mat(Kp));
    Dm = c_diff(i_mat(Km));
    % penalty term
    h_perp=Le; % temporary!
    pen = C_pen * (Dp/h_perp + Dm/h_perp) /2;
    
    % build the local edge gradient matrices
    % [[phi]],{{D.grad(b).ne}} 
    %      = (phi+ - phi-)(D+ grad(b+).ne + D- grad(b-).ne)/2
    %      =  (phi+) D+ grad(b+).ne/2 
    %       + (phi+) D- grad(b-).ne/2
    %       - (phi-) D+ grad(b+).ne/2
    %       - (phi-) D- grad(b-).ne/2
    % [[b]],{{D.grad(phi).ne}} 
    %      = (b+ - b-)(D+ grad(phi+).ne + D- grad(phi-).ne)/2
    %      =  (b+) D+ grad(phi+).ne/2 
    %       + (b+) D- grad(phi-).ne/2
    %       - (b-) D+ grad(phi+).ne/2
    %       - (b-) D- grad(phi-).ne/2
    % compute row vector: n' * G (ne is already stored as a 1x2 row vector)
    % -/-
    row_grad_m = ne * grad{Km}(:,:,indm);
    % edge matrix for this side
    col_b_m = zeros(nvm,1); col_b_m(1:2) = Le/2;
    aux = -Dm/2 * (col_b_m * row_grad_m + row_grad_m' * col_b_m');
    aux(1:2,1:2) = aux(1:2,1:2) + pen * Le * m1d;
    for i=1:nvm
        for j=1:nvm
            edgmat_mm(skip_m(i),skip_m(j)) = aux(i,j);
        end
    end
    A(gm(:),gm(:)) = A(gm(:),gm(:)) + edgmat_mm;

    % +/+
    row_grad_p = ne * grad{Kp}(:,:,indp);
    % edge matrix for this side
    col_b_p = zeros(nvp,1); col_b_p(1:2) = Le/2;
    aux = +Dp/2 * (col_b_p * row_grad_p + row_grad_p' * col_b_p');
    aux(1:2,1:2) = aux(1:2,1:2) + pen * Le * m1d;
    for i=1:nvp
        for j=1:nvp
            edgmat_pp(skip_p(i),skip_p(j)) = aux(i,j);
        end
    end
    A(gp(:),gp(:)) = A(gp(:),gp(:)) + edgmat_pp;

%     Dp=0;Dm=0;pen=-1;
%     A(:,:)=0;Le=6;
    pen=0;
    A(:,:)=0;
    % -(test)/+(solution)
    aux = ( -col_b_m * Dp/2*row_grad_p + Dm/2*row_grad_m' * col_b_p');
    aux(1:2,1:2) = aux(1:2,1:2) - pen * Le * m1d_mod;
    for i=1:nvm
        for j=1:nvp
            edgmat_mp(skip_m(i),skip_p(j)) = aux(i,j);
        end
    end
    A(gm(:),gp(:)) = A(gm(:),gp(:)) + edgmat_mp;

    % +(test)/-(solution)
    aux = ( col_b_p * Dm/2*row_grad_m + Dp/2*row_grad_p' * col_b_m');
    aux(1:2,1:2) = aux(1:2,1:2) - pen * Le * m1d_mod;
    for i=1:nvp
        for j=1:nvm
            edgmat_pm(skip_p(i),skip_m(j)) = aux(i,j);
        end
    end
    A(gp(:),gm(:)) = A(gp(:),gm(:)) + edgmat_pm;

    
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
