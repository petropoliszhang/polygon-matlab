function rectangular_pwld_driver()

%------------------------------------------------
close all; clc
% close all; clc; clear A; clear MM; clc
global verbose
verbose=false;
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

geofile='..\geom_codes\figs\random_quad_mesh_L100_n10_a0.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n10_a0.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n10_a0.txt';
%
geofile='..\geom_codes\figs\random_quad_mesh_L100_n10_a0.2.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n10_a0.2.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n10_a0.2.txt';
%
% geofile='..\geom_codes\figs\random_quad_mesh_L100_n10_a0.33.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_x_L100_n10_a0.33.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_mid_xy_L100_n10_a0.33.txt';
%
% geofile='..\geom_codes\figs\random_quad_mesh_L100_n30_a0.33.txt';
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

% geofile='..\geom_codes\figs\shestakov_quad_nc5_a0.25.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc4_a0.35.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc4_a0.25.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc4_a0.5.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc1_a0.25.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_L1_n2_a0.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc5_a0.25.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc6_a0.15.txt';
% geofile='..\geom_codes\figs\shestakov_quad_nc6_a0.25.txt';
% geofile='..\geom_codes\figs\random_quad_mesh_L100_n30_a0.33.txt';

% geofile='..\geom_codes\figs\random_poly_mesh_L10_n20_a0.95.txt';
% geofile='..\geom_codes\figs\random_poly_mesh_L10_n30_a0.95.txt';
%
% geofile='..\geom_codes\figs\smooth_poly_mesh_L10_n10_a0.8.txt';
% geofile='..\geom_codes\figs\smooth_poly_mesh_L10_n30_a0.8.txt';
% geofile='..\geom_codes\figs\smooth_poly_mesh_L10_n50_a0.88.txt';
% geofile='..\geom_codes\figs\smooth_poly_mesh_L1_n30_a0.1.txt';

% geofile='..\geom_codes\figs\random_quad_mesh_L1_n40_a0.txt';

geofile='..\geom_codes\figs\random_quad_mesh_L1_n2_a0.txt';

logi_mms  = true;
max_ref_cycles=15;
frac_ref=0.2;
mms_type=2;
logi_plot = true;
logi_plot_err_i = false;
generate_vtk_output = false;
vtk_basename = 'testing';
%
tot = 1/3; sca = 1/3;
c_diff=1/(3*tot); sigma_a=tot-sca; S_ext=0.10*0;
% bc type: 0= Dirichlet, homogeneous
%          1= Dirichlet, inhomogeneous
%          2= Neumann, homogeneous
%          3= Neumann, inhomogeneous
%          4= Robin phi/4 + D/2 \partial_n phi = Jinc
% values entered as LRBT
bc_type=[1 1 2 2  ];
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
C_pen=4*5;
C_pen_bd=2*C_pen;
%
%------------------------------------------------
%
% load mesh
%
fprintf('\n--------------------------------------------\n');
fprintf('------ Initial mesh (cycle 1) ------\n');
fprintf('--------------------------------------------\n');

[Lx,Ly,nel,ndof,connectivity,vert,n_edge,edg2poly,edg2vert,edg_perp,i_mat,i_src] =...
    read_geom(geofile);
% assign bc markers
edg2poly = assign_bc_markers(n_edge,edg2poly,edg2vert,vert,Lx,Ly);
% compute normal vectors
edg_normal = compute_edge_normals(n_edge,edg2vert,vert);
% assign the current refinement level
if(max_ref_cycles>1)
    curr_ref_lev = zeros(nel,1);
    corner_pos=cell(nel,1);
    for iel=1:nel
        corner_pos{iel}=1:4;
    end
end
%
%------------------------------------------------
%
% mms
%
if(logi_mms)
    bc_type=[0 0 0 0]; % imposed homogeneous Dirchlet
    switch(mms_type)
        case{0}
            mms=@(x,y)  S_ext+0*(x.*y);
            % select quadrature order
            n_quad = 8;
        case{1}
            % exact solution
            freq=1;
            exact=@(x,y) sin(freq*pi*x/Lx).*sin(freq*pi*y/Ly);
            % forcing rhs
            mms=@(x,y) (c_diff*(freq*pi)^2*(1/Lx^2+1/Ly^2)+sigma_a)*sin(freq*pi*x/Lx).*sin(freq*pi*y/Ly);
            % mms=@(x,y)  S_ext+0*(x.*y);
            % select quadrature order
            n_quad = 8;
        case{2}
            % exact solution
            x0=Lx*0.6;
            y0=Ly*0.7;
            varia=Lx^2/100;
            exact=@(x,y) 100*x.*(Lx-x).*y.*(Ly-y).*exp(-((x-x0).^2+(y-y0).^2)/varia)/(Lx*Ly)^2;
            % forcing rhs
            mms=@(x,y) c_diff.*(1.0./Lx.^2.*1.0./Ly.^2.*x.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(Lx-x).*2.0e2+1.0./Lx.^2.*1.0./Ly.^2.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(Ly-y).*2.0e2-(1.0./Lx.^2.*1.0./Ly.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(x.*2.0-x0.*2.0).*(Ly-y).*2.0e2)./varia-(1.0./Lx.^2.*1.0./Ly.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(y.*2.0-y0.*2.0).*(Lx-x).*2.0e2)./varia+(1.0./Lx.^2.*1.0./Ly.^2.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(x.*2.0-x0.*2.0).*(Lx-x).*(Ly-y).*2.0e2)./varia+(1.0./Lx.^2.*1.0./Ly.^2.*x.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(y.*2.0-y0.*2.0).*(Lx-x).*(Ly-y).*2.0e2)./varia+(1.0./Lx.^2.*1.0./Ly.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(Lx-x).*(Ly-y).*4.0e2)./varia-1.0./Lx.^2.*1.0./Ly.^2.*1.0./varia.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(x.*2.0-x0.*2.0).^2.*(Lx-x).*(Ly-y).*1.0e2-1.0./Lx.^2.*1.0./Ly.^2.*1.0./varia.^2.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(y.*2.0-y0.*2.0).^2.*(Lx-x).*(Ly-y).*1.0e2)+1.0./Lx.^2.*1.0./Ly.^2.*sigma_a.*x.*y.*exp(-((x-x0).^2+(y-y0).^2)./varia).*(Lx-x).*(Ly-y).*1.0e2;
            % select quadrature order
            n_quad = 8;
            xx=linspace(0,Lx);
            yy=linspace(0,Ly);
            [uu,vv]=meshgrid(xx,yy);
            zz=exact(uu,vv);
            figure(99)
            surf(uu,vv,zz)
            %             figure(999)
            %             surf(uu,vv,mms(uu,vv))
        otherwise
            error('wrong mss type');
    end
else
    exact='';
    mms='';
    n_quad=0;
end

%=============================================
%=============================================
% refinement cycles
%=============================================
%=============================================
for i_cycle=1:max_ref_cycles
    
    if(i_cycle==1)
        t_cycle_beg = t_beg;
    else
        t_cycle_beg = cputime;
    end
    %------------------------------------------------
    %
    % assemble + solve
    %
    z = DG_assemble_solve( ndof,nel,n_edge,vert,connectivity,edg2poly,edg2vert,edg_normal,edg_perp,C_pen,C_pen_bd,...
        i_mat,c_diff,sigma_a,i_src,S_ext,logi_mms,mms,n_quad,bc_type,bc_val );
    
    %------------------------------------------------
    %
    % error indicator
    %
    err_i = error_ind(z,nel,n_edge,vert,connectivity,edg2poly,edg2vert,c_diff,i_mat);
    
    if(logi_plot_err_i)
        incr=0;
        figure(13+(i_cycle-1)*10);clf
        for iel=1:nel
            g=connectivity{iel}(:);
            ee=log10(err_i(iel)*ones(length(g),1));
            patch(vert(g,1),vert(g,2),ee,ee,'FaceColor','interp'); %,'LineStyle','none');
        end
        view(-135,25);
        view(0,90);
        figure(14+(i_cycle-1)*10);clf
        for iel=1:nel
            g=connectivity{iel}(:);
            ee=(err_i(iel)*ones(length(g),1));
            patch(vert(g,1),vert(g,2),ee,ee,'FaceColor','interp'); %,'LineStyle','none');
        end
        view(-135,25);
        view(0,90);
    end

    %------------------------------------------------
    %
    % plot
    %
    if(logi_plot)
        
        incr=0;
        figure(11+(i_cycle-1)*incr);clf
        for iel=1:nel
            g=connectivity{iel}(:);
            patch(vert(g,1),vert(g,2),z(g),z(g),'FaceColor','interp'); %,'LineStyle','none');
        end
        view(-135,25);
        view(0,90);
        figure(12+(i_cycle-1)*incr);clf
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
        
        norm_data(i_cycle,1)=ndof;
        norm_data(i_cycle,2)=L2_norm(ndof,nel,connectivity,vert,n_quad,z,exact);
        
    end % end logical test
    
    %------------------------------------------------
    %
    % vtk output
    %
    if generate_vtk_output
        cycle_number = i_cycle;
        if(max_ref_cycles==1),  cycle_number=[]; end
        create_vtk_output(vtk_basename,ndof,nel,connectivity,vert,z,cycle_number);
    end
    
    %------------------------------------------------
    t_cycle_end = cputime;
    fprintf('\n\nTime in cycle %d = %g \n',i_cycle,t_cycle_end-t_cycle_beg);
    %------------------------------------------------
    %
    % refinement
    if(i_cycle<max_ref_cycles)
        
        fprintf('\n--------------------------------------------\n');
        fprintf('------ Refinement cycle # %3.3d ------\n',i_cycle+1);
        fprintf('--------------------------------------------\n');
        
        %         % for debugging only !
        %         err_i(:)=0;
        % %         err_i(end)=1;
        %         distance_=zeros(nel,1);
        %         for iel=1:nel
        %             g = connectivity{iel}(corner_pos{iel});
        %             v = vert(g,:);
        %             vc=mean(v);
        %             distance_(iel)=norm(vc-[Lx Ly]);
        %         end
        %         [dummy,ind_]=min(distance_);
        %         err_i(ind_)=1;
        
%         figure(99);clf
%         hold all
%         for id=1:length(connectivity)
%             vv= vert(connectivity{id},:);
%             xc=mean(vv(:,1));
%             yc=mean(vv(:,2));
%             vv(end+1,:)=vv(1,:);
%             plot(vv(:,1),vv(:,2),'+-')
%             %plot(xc,yc,'x');
%             str=sprintf('%d',id);
%             text(xc,yc,str);
%             %     pause
%         end
        
        [nel,ndof,connectivity,corner_pos,vert,n_edge,edg2poly,edg2vert,edg_perp,i_mat,i_src,curr_ref_lev] =...
            refine_geom_v2(err_i,frac_ref,curr_ref_lev,nel,connectivity,corner_pos,vert,...
            edg2poly,edg2vert,i_mat,i_src);
        
        % assign bc markers
        edg2poly = assign_bc_markers(n_edge,edg2poly,edg2vert,vert,Lx,Ly);
        
        % compute normal vectors
        edg_normal = compute_edge_normals(n_edge,edg2vert,vert);
        
    end
    
    %------------------------------------------------
    
end
%=============================================
%=============================================
% refinement cycles
%=============================================
%=============================================
t_end=cputime;
fprintf('\n\n-----------------------------\nTotal time    = %g \n',t_end-t_beg);

if(logi_mms)
    figure(999)
    plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-')
end

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


