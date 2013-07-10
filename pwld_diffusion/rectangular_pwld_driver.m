function rectangular_pwld()

%------------------------------------------------
close all; clc
% close all; clc; clear A; clear MM; clc
%------------------------------------------------
%
% data
%
geofile='..\geom_codes\figs\random_quad_mesh_L100_n50_a0.33.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n50_a0.33.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n50_a0.33.txt';

% geofile='..\geom_codes\figs\random_quad_mesh_L100_n50_a0.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n50_a0.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n50_a0.txt';

% geofile='..\geom_codes\figs\random_quad_mesh_L100_n10_a0.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n10_a0.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n10_a0.txt';
% 
% geofile='..\geom_codes\figs\random_quad_mesh_L100_n10_a0.2.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n10_a0.2.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n10_a0.2.txt';
% 
% geofile='..\geom_codes\figs\random_quad_mesh_L100_n10_a0.33.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n10_a0.33.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n10_a0.33.txt';
% 
geofile='..\geom_codes\figs\random_quad_mesh_L100_n30_a0.33.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n30_a0.33.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n30_a0.33.txt';

% geofile='..\geom_codes\figs\random_quad_mesh_L100_n30_a0.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n30_a0.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n30_a0.txt';
% 
% geofile='..\geom_codes\figs\random_quad_mesh_L100_n30_a0.1.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n30_a0.1.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n30_a0.1.txt';

% geofile='..\geom_codes\figs\random_quad_mesh_L100_n2_a0.1.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n2_a0.1.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n2_a0.1.txt';

% geofile='..\geom_codes\figs\random_quad_mesh_L100_n3_a0.25.txt';

geofile='..\geom_codes\figs\shestakov_quad_nc5_a0.25.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc4_a0.35.txt';
geofile='..\geom_codes\figs\shestakov_quad_nc4_a0.25.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc4_a0.5.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc1_a0.25.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_L1_n2_a0.txt';


geofile='..\geom_codes\\figs\smooth_quad_mesh_L100_n30_a0.15.txt';
geofile='..\geom_codes\\figs\smooth_quad_mesh_L100_n50_a0.15.txt';

%
logi_mms  = true;
logi_plot = true;
vtk_basename = 'rectangular';
%
tot = 1/3; sca = 1/3;
c_diff=1/(3*tot); sigma_a=tot-sca; S_ext=1.0; 
% bc type: 0= Dirichlet, homogeneous
%          1= Dirichlet, inhomogeneous
%          2= Neumann, homogeneous
%          3= Neumann, inhomogeneous
%          4= Robin phi/4 + D/2 \partial_n phi = Jinc
% values entered as LRBT
bc_type=[0 0 0 0 ];
bc_val.left  = 100;
bc_val.right = -50;
bc_val.bottom= 50;
bc_val.top   = 10;
%
%------------------------------------------------
t_beg=cputime;
%------------------------------------------------
%
% numerical parameters
%
C_pen=4;
C_pen_bd=2*C_pen;
%
%------------------------------------------------
%
% load mesh 
%
[Lx,Ly,nel,ndof,connectivity,vert,n_edge,edg2poly,edg2vert,i_mat,i_src] =...
    read_geom(geofile);
% assign bc markers
edg2poly = assign_bc_markers(n_edge,edg2poly,edg2vert,vert,Lx,Ly);
% compute normal vectors
edg_normal = compute_edge_normals(n_edge,edg2vert,vert);
%
%------------------------------------------------
%
% mms
% 
if(logi_mms)
    bc_type=[0 0 0 0]; % imposed homogeneous Dirchlet
    % exact solution
    freq=3;
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
% assemble + solve
%
z = DG_assemble_solve( ndof,nel,n_edge,vert,connectivity,edg2poly,edg2vert,edg_normal,C_pen,C_pen_bd,...
    i_mat,c_diff,sigma_a,i_src,S_ext,logi_mms,mms,n_quad,bc_type,bc_val );

%------------------------------------------------
%
% plot
%
if(logi_plot)
    
    figure(11);clf
    for iel=1:nel
        g=connectivity{iel}(:);
        patch(vert(g,1),vert(g,2),z(g),z(g),'FaceColor','interp'); %,'LineStyle','none');
    end
    view(-135,25);
    figure(12);clf
    % plot on finer mesh
    % 4---3   vertex anti-clockwise ordering,
    % | c |
    % 1---2
    for iel=1:nel
        g=connectivity{iel}(:);
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
%
% L-2 norm
%
if(logi_mms)

    L2_norm(ndof,nel,connectivity,vert,n_quad,z,exact);

end % end logical test

%------------------------------------------------
%
% vtk output 
%
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


