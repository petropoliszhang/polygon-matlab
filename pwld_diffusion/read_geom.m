function [Lx,Ly,nel,ndof,connectivity,vert,n_edge,edg2poly,edg2vert,edg_perp,i_mat,i_src] ...
    = read_geom(file)

t1=cputime;

fid=fopen(file);
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

% complete mesh data
[n_edge,edg2poly,edg2vert,edg_perp ] = complete_mesh_data( nel,ndof,n_grid_vert,vert,vert_grid,connectivity );

t2=cputime;
fprintf('Mesh time     = %g \n\n',t2-t1);

return
end
