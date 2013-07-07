function rectangular_pwld()

%------------------------------------------------
close all; clc; clear A; clear MM;
%------------------------------------------------
t_beg=cputime;
%------------------------------------------------
% clear all; close all; clc
%
% data
%
tot=1/3;sca=1/3;
Lx=100; c_diff=1/(3*tot); sigma_a=tot-sca; S_ext=0.10; Ly=Lx;
% bc type: 0= Dirichlet, homogeneous
%          1= Dirichlet, inhomogeneous
%          2= Neumann, homogeneous
%          3= Neumann, inhomogeneous
%          4= Robin phi/4 + D/2 \partial_n phi = Jinc
% values entered as LRBT
bc_type=[0 0 0 0];
bc_val.left  = 0;
bc_val.right = 0;
bc_val.bottom= 10;
bc_val.top   = 0;
%
logi_plot = true;
vtk_basename = 'rectangular';
%
% numerical parameters
%
C_pen=4;
C_pen_bd=1*C_pen;
%
%------------------------------------------------
logi_mms=true;
if(logi_mms)
    bc_type=[0 0 0 0]; % imposed homogeneous Dirchlet
    % exact solution
    freq=1;
    exact=@(x,y) sin(freq*pi*x/Lx).*sin(freq*pi*y/Ly);
    % forcing rhs
    mms=@(x,y) (c_diff*(freq*pi)^2*(1/Lx^2+1/Ly^2)+sigma_a)*sin(freq*pi*x/Lx).*sin(freq*pi*y/Ly);
    % mms=@(x,y)  S_ext+0*(x.*y);
    % select quadrature order
    n_quad = 8;
else
    exact='';
    mms='';
    n_quad=0;
end
%------------------------------------------------
%
% load mesh 
%
nx=2^4; ny=nx;
x=linspace(0,Lx,nx+1); y=linspace(0,Ly,ny+1);
nel=nx*ny;
ndof = 4*nel;
i_mat=ones(nel,1);
%
%------------------------------------------------
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
%------------------------------------------------
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

%------------------------------------------------
edg2poly = assign_bc_markers(n_edge,edg2poly,edg2vert,vert,Lx,Ly);
%------------------------------------------------
edg_normal = compute_edge_normals(n_edge,edg2vert,vert);

%------------------------------------------------
% assemble + solve
z = DG_assemble_solve( ndof,nel,n_edge,vert,connectivity,edg2poly,edg2vert,edg_normal,C_pen,C_pen_bd,...
    i_mat,c_diff,sigma_a,S_ext,logi_mms,mms,n_quad,bc_type,bc_val );

%------------------------------------------------
% plot
if(logi_plot)
    
    figure(11);clf
    for iel=1:nel
        g=connectivity(iel,:);
        patch(vert(g,1),vert(g,2),z(g),z(g),'FaceColor','interp'); %,'LineStyle','none');
    end
    view(-135,25);
    figure(12);clf
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
    view(-135,25);

end

%------------------------------------------------
% L-2 norm
if(logi_mms)

    L2_norm(ndof,nel,connectivity,vert,n_quad,z,exact);

end % end logical test

%------------------------------------------------
% vtk output 
create_vtk_output(vtk_basename,ndof,nel,connectivity,vert,z)

%------------------------------------------------
t_end=cputime;
fprintf('\n\n-----------------------------\nTotal time    = %g \n',t_end-t_beg);

return
end
%------------------------------------------------


%------------------------------------------------
% linear 1d solution with robin on the left and right with no volumetric
% source and absorption=0
%
% phi(x) = a.x + b
%
% bc left : phi/4(0) - D/2 dphi/dx|_0 = J
% bc right: phi/4(L) + D/2 dphi/dx|_L = 0
%
% a = -4J/(L+4D)
% b =  4J(L+2D)/(L+4D)
%
% phi(0) = b
% phi(L) = a.L + b = 8JD/(L+4D)
%


