% compute the gradient per side of an arbitrary polygon
clear all; close all; clc;

tri=false;
if(tri)
    poly=[1 2 3]
    vert=[0 0; 1 0; 0 1]
else
    poly=[1 2 3 4]
    vert=[0 0; 1 0; 1 1; 0 1]
end

% centroid
C=mean(vert);
% alpha coef
nv=length(poly);
alpha=1/nv;

g=zeros(2,nv,nv);
% list of vertices, when looping over sides
list_vert=1:nv;

% loop over sides
for iside=1:nv
    % pick 1st vertex
    irow1=iside;
    % pick 2ndt vertex
    irow2=irow1+1; if(irow2>nv), irow2=1; end
    % compute J=2.Area for that side
    A=vert(irow1,:); B=vert(irow2,:);
    Jac=[(B-A)' (C-A)'];
    iJt=inv(Jac');
    % stiffness matrix
    g_side=[-1 1 0; (alpha-1) alpha alpha];
    aux=iJt*g_side;
    g(:,list_vert(1:2),iside) = aux(1:2,1:2);
    for k=3:nv
        g(:,list_vert(k),iside) = aux(1:2,3);
    end
    % shift vertex IDs
    list_vert=[list_vert list_vert(1)];
    list_vert(1)=[];
end
g
