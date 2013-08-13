clear all; close all; clc

load .\results\convergence\poly\shes\shes_poly_mms_1_L1_nc6_a0.1.mat
n_vertices=zeros(nel,1);

for iel=1:nel
    n_vertices(iel)=length(connectivity{iel});
end


min_vert=min(n_vertices);
max_vert=max(n_vertices);
for k=min_vert:max_vert
    nnn=find(n_vertices==k);
    fprintf('Number of polygons with %d vertices = %d \n',k,length(nnn));
end
fprintf('Total number of polygons           = %d \n',nel);

nnn=find(n_vertices==max_vert)

figure(11)
hold all
for id=1:length(connectivity)
    ind=find(nnn==id);
    if isempty(ind)
        patch(vert(connectivity{id},1),vert(connectivity{id},2),'white'); % use color i.
    else
        disp('here')
        patch(vert(connectivity{id},1),vert(connectivity{id},2),'red'); % use color i.
    end
end
set(gca,'PlotBoxAspectRatio',[5 5 1])


clear all; 

load .\results\convergence\poly\z-mesh\z_mesh_poly_n80_a0.05_mms_1.mat
n_vertices=zeros(nel,1);

for iel=1:nel
    n_vertices(iel)=length(connectivity{iel});
end


min_vert=min(n_vertices);
max_vert=max(n_vertices);
for k=min_vert:max_vert
    nnn=find(n_vertices==k);
    fprintf('Number of polygons with %d vertices = %d \n',k,length(nnn));
end
fprintf('Total number of polygons           = %d \n',nel);

nnn=find(n_vertices==max_vert)

figure(12)
hold all
for id=1:length(connectivity)
    ind=find(nnn==id);
    if isempty(ind)
        patch(vert(connectivity{id},1),vert(connectivity{id},2),'white'); % use color i.
    else
        disp('here')
        patch(vert(connectivity{id},1),vert(connectivity{id},2),'red'); % use color i.
    end
end
set(gca,'PlotBoxAspectRatio',[5 5 1])
