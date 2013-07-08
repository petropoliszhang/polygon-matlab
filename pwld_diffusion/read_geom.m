function [Lx,Ly,nel,ndof,connectivity,vert,n_edge,edg2poly,edg2vert,i_mat,i_src] = read_geom(file);

t1=cputime;

% clear all; close all; clc;
% file='.\geom_codes\figs\aa_random_quad_mesh_L1_n10_a0.53.txt'

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

% read DG vertices (counter-clockwise)
n_vert=C(ind);
ind=ind+1;
for i=1:n_vert
    vert(i,1)=C(ind);
    ind=ind+1;
    vert(i,2)=C(ind);
    ind=ind+1;
end

% read grid vertices (counter-clockwise)
n_grid_vert=C(ind);
ind=ind+1;
for i=1:n_grid_vert
    vert_grid(i,1)=C(ind);
    ind=ind+1;
    vert_grid(i,2)=C(ind);
    ind=ind+1;
end

% close geom data file
fclose(fid);

% find relationship between vert and vert_grid
vert_link=zeros(ndof,1);
for i=1:ndof
    v=vert(i,1:2);
    % get difference
    aux = vert_grid - kron(ones(length(vert_grid(:,1)),1),v);
    % get vector of norm
    aa=sqrt(aux(:,1).^2+aux(:,2).^2);
    ind=find(aa<1e-12);
    len=length(ind);
    if(len==0 || len >1), 
        error('vert_link'); 
    end
    vert_link(i)=ind(1);
end

% edge data
n_edge=0;
edg2poly=zeros(0,2);
edg_normal=zeros(0,2);
edg2vert=zeros(0,4);
for iel=1:nel
    elem=connectivity{iel}(:);
    nedg=length(elem);
    elem(end+1)=elem(1);
    for i=1:nedg
        ed=elem(i:i+1);
        [edg2poly,edg2vert,n_edge]=is_edge_already_recorded(...
            ed,edg2poly,edg2vert,iel,vert_link,n_edge);
    end
end
clear vert_link; % not needed any longer

% check orientation, centroid
tot_area=0;
for iel=1:nel
    g=connectivity{iel}(:);
    xx=vert(g,1); yy=vert(g,2);
    
    % check orientation, verify area
    [or,ar] = polyorient(xx,yy);
    if(or~=1), warning('orientation problem'); end
    tot_area=tot_area+ar;
    
    % check centroid
    xc=mean(xx); yc=mean(yy);
    [in on]=inpolygon(xc,yc,xx,yy);
    % find in points
    ind_in=find(in==1);
    % find if these 'in' points are on the poly edge
    in_and_on = on(ind_in);
    ind_in_and_on=find(in_and_on==1);
    if(~isempty(ind_in_and_on))
        warning('centroid is on poly edge');
    end
    % find out points
    ind_out=find(in==0);
    if(~isempty(ind_out))
        error('centroid is OUTSIDE of poly');
    end
end
fprintf('tot_area = %g \n',tot_area);

t2=cputime;
fprintf('Mesh time     = %g \n',t2-t1);

return
end