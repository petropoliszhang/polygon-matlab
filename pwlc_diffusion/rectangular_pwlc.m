clear all; close all; clc
% data
Lx=1; D=1; S=0; Q=0; Ly=1
% geometry
nx=10; ny=nx;
x=linspace(0,Lx,nx+1); y=linspace(0,Ly,ny+1);
nel=nx*ny;
ndof = (nx+1)*(ny+1);
% 4---3   vertex anti-clockwise ordering,
% |   |
% 1---2
connectivity=zeros(nel,4);
iel=0;
for j=1:ny
    for i=1:nx
        iel=iel+1;
        i1=(j-1)*(nx+1)+i; i2=i1+1;
        i4=(j  )*(nx+1)+i; i3=i4+1;
        connectivity(iel,:)=[i1 i2 i3 i4];
    end
end
% coordinates
ind=0;
vert=zeros(ndof,2);
for j=1:ny+1
    for i=1:nx+1
        ind=ind+1;
        vert(ind,1:2)=[x(i) y(j)];
    end
end

% assemble
A = spalloc(ndof,ndof,9); b=zeros(ndof,1);
for iel=1:nel
    g=connectivity(iel,:);
    v=vert(g,:);
    [M,K,f]=build_pwld_local_matrices(g,v);
    A(g(:),g(:)) = A(g(:),g(:)) + D*K +S*M;
    b(g(:)) = b(g(:)) + Q*f;
end
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

