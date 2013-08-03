function [Lx,Ly,nel,ndof,connectivity,vert,n_edge,edg2poly,edg2vert,edg_perp,i_mat,i_src] ...
    = read_geom(file)

t1=cputime;

%------------------------------------------------
fid=fopen(file);

%------------------------------------------------
if(fid<0)
    %------------------------------------------------
    % the file does not exist, we create a uniform mesh directly
    fprintf('The Geo file does not exist. A uniform mesh will be created instead\n');
    Lx=1;Ly=Lx;
    nx=input('Enter the # of subdivisions along x (a 0 or negative value will abort): ');
    fprintf('\n');
    if(nx<=0)
        error('Aborting. nx<=0');
    end
    ny=nx;
    t1=cputime; % reset
    x=linspace(0,Lx,nx+1); y=linspace(0,Ly,ny+1);
    nel=nx*ny;
    i_mat=ones(nel,1);
    i_src=ones(nel,1);
    ndof = 4*nel;
    
    %------------------------------------------------
    % 4---3   vertex anti-clockwise ordering,
    % |   |
    % 1---2
    connectivity=cell(nel,1);
    for iel=1:nel
        skip = 4*(iel-1);
        i1 = skip + 1;
        i2 = skip + 2;
        i3 = skip + 3;
        i4 = skip + 4;
        connectivity{iel} = zeros(4,1);
        connectivity{iel}(:)=[i1 i2 i3 i4];
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
    n_grid_vert = length( vert_grid(:,1) );
%     % find relationship between vert and vert_grid
%     vert_link=zeros(ndof,1);
%     for i=1:ndof
%         v=vert(i,1:2);
%         for k=1:(nx+1)*(ny+1)
%             if(norm( vert_grid(k,:)-v )<1e-12)
%                 if(vert_link(i)~=0), error('vert_link(i) should be 0'); end
%                 vert_link(i)=k;
%             end
%         end
%     end
    
    %------------------------------------------------

else
    
    %%% ------------------  READ FILE
    a=textscan(fid,'%n','commentstyle','#');
    C=a{1};clear a;
    
    % get dimensions
    Lx = C(1);
    Ly = C(2);
    
    % get # of polygons
    nel = C(3);
    
    % initialize arrays
    connectivity=cell(nel,1);
    n_vertices=zeros(nel,1);
    i_mat=zeros(nel,1);
    i_src=zeros(nel,1);
    
    % read connectivity
    ind=4;
    for iel=1:nel
        n_vertices(iel)=C(ind);
        i1=ind+1;
        i2=ind+n_vertices(iel);
        connectivity{iel} = zeros(n_vertices(iel),1);
        connectivity{iel}(:)=C(i1:i2);
        ind=i2+1;
        i_mat(iel)=C(ind);
        ind=ind+1;
        i_src(iel)=C(ind);
        ind=ind+1;
    end
    
    ndof=sum(n_vertices);
    
    min_vert=min(n_vertices);
    max_vert=max(n_vertices);
    for k=min_vert:max_vert
        nnn=find(n_vertices==k);
        fprintf('Number of polygons with %d vertices = %d \n',k,length(nnn));
    end
    fprintf('Total number of polygons           = %d \n',nel);
    
    % read DG vertices (counter-clockwise)
    n_vert=C(ind);
    vert=zeros(n_vert,2);
    ind=ind+1;
    for i=1:n_vert
        vert(i,1)=C(ind);
        ind=ind+1;
        vert(i,2)=C(ind);
        ind=ind+1;
    end
    
    % read grid vertices
    n_grid_vert=C(ind);
    vert_grid=zeros(n_grid_vert,2);
    ind=ind+1;
    for i=1:n_grid_vert
        vert_grid(i,1)=C(ind);
        ind=ind+1;
        vert_grid(i,2)=C(ind);
        ind=ind+1;
    end
    
    % close geom data file
    fclose(fid);
end
%%% ------------------ END READ file
%------------------------------------------------

% convexity
convexity_util(connectivity,nel,vert);

% complete mesh data
[n_edge,edg2poly,edg2vert,edg_perp ] = complete_mesh_data( nel,ndof,n_grid_vert,vert,vert_grid,connectivity );

t2=cputime;
fprintf('Mesh time     = %g \n\n',t2-t1);

return
end
